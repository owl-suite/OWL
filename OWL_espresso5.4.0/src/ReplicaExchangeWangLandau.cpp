#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>             // std::string
#include <sstream>            // std::istringstream
#include "ReplicaExchangeWangLandau.hpp"
#include "RandomNumberGenerator.hpp"

// Constructor
ReplicaExchangeWangLandau::ReplicaExchangeWangLandau(MPICommunicator PhySystemComm, MPICommunicator MCAlgorithmComm) : h(simInfo.restartFlag, simInfo.MCInputFile)
{

  std::cout << "Simulation method: Replica-Exchange Wang-Landau sampling\n";

  /// Pass MPI communicators from arguments
  PhysicalSystemComm = PhySystemComm;
  REWLComm = MCAlgorithmComm; 

  /// Set initial values for private members
  numWalkers = simInfo.numWalkers;
  readREWLInputFile(simInfo.MCInputFile);

  /// group processors into different walkers
  walkerID = (GlobalComm.thisMPIrank - (GlobalComm.thisMPIrank % simInfo.numMPIranksPerWalker)) / simInfo.numMPIranksPerWalker;
 
  myWindow = ( walkerID - (walkerID % numWalkersPerWindow) ) / numWalkersPerWindow;

  //YingWai's check
  //printf("YingWai's check: Inside REWL constructor. numWalkers = %3d, world_rank = %3d, myWindow = %3d, walkerID = %3d, numWindows = %3d\n", numWalkers, GlobalComm.thisMPIrank, myWindow, walkerID, numWindows);

  partnerID = -1;
  partnerWindow = -1;

  swapDirection = 0;

  upExchanges = 0;
  downExchanges = 0;

  GlobalComm.barrier();

}


ReplicaExchangeWangLandau::~ReplicaExchangeWangLandau()
{

  //printf("Exiting ReplicaExchangeWangLandau class... \n");
  //delete physical_system;

}


// Private member functions
void ReplicaExchangeWangLandau::assignSwapPartner()
{

  partnerID = -1; 
  partnerWindow = -1; 

  // Find the partner's window
  switch (swapDirection) {
  case 0 :                              // 0-1 & 2-3 & .... & [numprocs-1]-nobody
  {
    if ((myWindow % 2 == 0) && (myWindow < (numWindows - 1))) {
      partnerWindow = myWindow + 1;
      upExchanges++;
    }   
    else if (myWindow % 2 != 0) {
      partnerWindow = myWindow - 1;
      downExchanges++;
    }   
    swapDirection = 1;
    //printf("YingWai's check: Walker. %5d, partner is %5d, swapDirection is %5d\n", myWindow, partnerWindow, swapDirection);
    break;
  }
  case 1 :                              // 0-nobody & 1-2 & 3-4 & ....
  {
    if ((myWindow % 2 == 0) && (myWindow != 0)) {
      partnerWindow = myWindow - 1;
      downExchanges++;
    }   
    else if ((myWindow % 2 != 0) && (myWindow < (numWindows - 1))) {
      partnerWindow = myWindow + 1;
      upExchanges++;
    }   
    swapDirection = 0;
    //printf("YingWai's check: Walker. %5d, partner is %5d, swapDirection is %5d\n", myWindow, partnerWindow, swapDirection);
  }
  }

  // Find the partner's ID
  if (numWalkersPerWindow == 1) {
    partnerID = partnerWindow;
  }
  else {  // multiple walkers per window case
    // to be implemented...
  }

}


void ReplicaExchangeWangLandau::swapEnergy(double &energyForSwap)
{

  if (partnerID != -1)       // Swap energy with partner  
    REWLComm.swapScalar(energyForSwap, partnerID);

}


bool ReplicaExchangeWangLandau::determineAcceptance(double myDOSRatio)
{

  double partnerDOSRatio {0.0};
  int change {0};

  if (myWindow % 2 == 1) {          // Receiver and calculator 
    REWLComm.recvScalar(partnerDOSRatio, partnerID);

    double acceptProb = myDOSRatio * partnerDOSRatio;

    if (getRandomNumber2() < acceptProb) change = 1;

    REWLComm.sendScalar(change, partnerID);
  }
  else {                            // Sender and non-calculator
    REWLComm.sendScalar(myDOSRatio, partnerID);
    REWLComm.recvScalar(change, partnerID);
  }

  if (change) return true;
  else return false;

}


void ReplicaExchangeWangLandau::swapConfiguration(double evecsForSwap[], int numElements)
{
  if (partnerID != -1)
    REWLComm.swapVector(evecsForSwap, numElements, partnerID);

}


// Public member functions

void ReplicaExchangeWangLandau::run(PhysicalSystem* physical_system)
{

  //double currentTime, lastBackUpTime;
  currentTime = lastBackUpTime = MPI_Wtime();
  if (GlobalComm.thisMPIrank == 0)
    printf("Running ReplicaExchangeWangLandau...\n");

  char fileName[51];

  // Find the first energy that falls within the WL energy range    
  while (!acceptMove) {
    physical_system -> doMCMove();
    physical_system -> getObservables();
    physical_system -> acceptMCMove();    // always accept the move to push the state forward
    acceptMove = h.checkEnergyInRange(physical_system -> observables[0]);
  }

  // Always accept the first energy if it is within range
  h.updateHistogramDOS(physical_system -> observables[0]);

  // Write out the energy
  if (GlobalComm.thisMPIrank == 0) 
    physical_system -> writeConfiguration(0, "energyLatticePos.dat");

//-------------- End initialization --------------//

// WL procedure starts here
  while (h.modFactor > h.modFactorFinal) {
    h.histogramFlat = false;

    while (!(h.histogramFlat)) {
      for (unsigned int MCSteps=0; MCSteps<h.histogramCheckInterval; MCSteps++) {

        physical_system -> doMCMove();
        physical_system -> getObservables();

        // check if the energy falls within the energy range
        if ( !h.checkEnergyInRange(physical_system -> observables[0]) )
          acceptMove = false;
        else {
          // determine WL acceptance
          if ( exp(h.getDOS(physical_system -> oldObservables[0]) - 
                   h.getDOS(physical_system -> observables[0])) > getRandomNumber2() )
            acceptMove = true;
          else
            acceptMove = false;
        }

        if (acceptMove) {
           // Update histogram and DOS with trialEnergy
           h.updateHistogramDOS(physical_system -> observables[0]);
           h.acceptedMoves++;        
 
           physical_system -> acceptMCMove();
        }
        else {
           physical_system -> rejectMCMove();

           // Update histogram and DOS with oldEnergy
           h.updateHistogramDOS(physical_system -> oldObservables[0]);
           h.rejectedMoves++;
        }
        h.totalMCsteps++;

        // Write restart files at interval
        currentTime = MPI_Wtime();
        if (PhysicalSystemComm.thisMPIrank == 0) {
          if (currentTime - lastBackUpTime > 300) {
            sprintf(fileName, "hist_dos_checkpoint_walker%05d.dat", REWLComm.thisMPIrank);
            h.writeHistogramDOSFile(fileName, h.iterations, REWLComm.thisMPIrank);
            sprintf(fileName, "OWL_checkpoint_walker%05d.dat", REWLComm.thisMPIrank);
            physical_system -> writeConfiguration(1, fileName);
            lastBackUpTime = currentTime;
          }
        }
      }
      //h.writeHistogramDOSFile("hist_dos_checkpoint.dat");

      // Check histogram flatness
      h.histogramFlat = h.checkHistogramFlatness();
      //if (h.numHistogramNotImproved >= h.histogramRefreshInterval)
      //  h.refreshHistogram();
    }

    if (GlobalComm.thisMPIrank == 0)
      printf("Number of iterations performed = %d\n", h.iterations);
    
      // Also write restart file here 
    if (PhysicalSystemComm.thisMPIrank == 0) {
      sprintf(fileName, "hist_dos_iteration%02d_walker%05d.dat", h.iterations, REWLComm.thisMPIrank);
      h.writeHistogramDOSFile(fileName, h.iterations, REWLComm.thisMPIrank);
      sprintf(fileName, "OWL_checkpoint_walker%05d.dat", REWLComm.thisMPIrank);
      physical_system -> writeConfiguration(1, fileName);
    }

    // Go to next iteration
    h.modFactor /= h.modFactorReducer;
    h.resetHistogram();
    h.iterations++;
  }

  // Write out data at the end of the simulation
  sprintf(fileName, "dos_walker%05d.dat", REWLComm.thisMPIrank);
  h.writeNormDOSFile(fileName, REWLComm.thisMPIrank);
  sprintf(fileName, "hist_dos_final_walker%05d.dat", REWLComm.thisMPIrank);
  h.writeHistogramDOSFile(fileName, h.iterations, REWLComm.thisMPIrank);

}


void ReplicaExchangeWangLandau::readREWLInputFile(const char* fileName)
{
  
  std::cout << "Reading REWL input file: " << fileName << std::endl;

  std::ifstream inputFile(fileName);
  std::string line, key;

  if (inputFile.is_open()) {
    while (std::getline(inputFile, line)) {
  
      if (!line.empty()) {
  
        std::istringstream lineStream(line);
        lineStream >> key;
  
        if (key.compare(0, 1, "#") != 0) {
  
          // Read and set numWindows, numWalkersPerWindow, overlap
          if (key == "numberOfWindows") {
            lineStream >> numWindows;
            //std::cout << "REWL: numberOfWindows = " << numWindows << std::endl;
            continue;
          }
          if (key == "numberOfWalkersPerWindow") {
            lineStream >> numWalkersPerWindow;
            //std::cout << "REWL: numberOfWalkersPerWindow = " << numWalkersPerWindow << std::endl;
            continue;
          }
          if (key == "overlap") {
            lineStream >> overlap;
            //std::cout << "REWL: overlap = " << overlap << std::endl;
            continue;
          }
  
        }

      }

    }
    inputFile.close();
  }


}

