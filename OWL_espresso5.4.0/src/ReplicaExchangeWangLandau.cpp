#include <cstdio>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>             // std::string
#include <sstream>            // std::istringstream
#include "ReplicaExchangeWangLandau.hpp"
#include "RandomNumberGenerator.hpp"

// Constructor
ReplicaExchangeWangLandau::ReplicaExchangeWangLandau(PhysicalSystem* ps, MPICommunicator PhySystemComm, MPICommunicator MCAlgorithmComm) : h(simInfo.restartFlag, simInfo.MCInputFile)
{

  std::cout << "Simulation method: Replica-Exchange Wang-Landau sampling\n";

  physical_system = ps;

  /// Pass MPI communicators from arguments
  PhysicalSystemComm = PhySystemComm;
  REWLComm = MCAlgorithmComm; 

  //if (REWLComm.communicator != MPI_COMM_NULL)
  //  std::cout << "I am Global ID = " << GlobalComm.thisMPIrank << ", REWLComm ID = " << REWLComm.thisMPIrank << std::endl;

  /// Set initial values for private members
  numWalkers = simInfo.numWalkers;
  readREWLInputFile(simInfo.MCInputFile);

  /// group processors into different walkers
  walkerID = (GlobalComm.thisMPIrank - (GlobalComm.thisMPIrank % simInfo.numMPIranksPerWalker)) / simInfo.numMPIranksPerWalker;
 
  myWindow = ( walkerID - (walkerID % numWalkersPerWindow) ) / numWalkersPerWindow;

  //YingWai's check
  //printf("YingWai's check: Inside REWL constructor. numWalkers = %3d, world_rank = %3d, myWindow = %3d, walkerID = %3d, numWindows = %3d\n", numWalkers, GlobalComm.thisMPIrank, myWindow, walkerID, numWindows);

  partnerID     = -1;
  partnerWindow = -1;

  swapDirection = 0;

  upExchanges   = 0;
  downExchanges = 0;

  MaxModFactor  = std::numeric_limits<double>::max();

  GlobalComm.barrier();

}


ReplicaExchangeWangLandau::~ReplicaExchangeWangLandau()
{

  if (GlobalComm.thisMPIrank == 0)
    printf("Exiting ReplicaExchangeWangLandau class... \n");

}

/////////////////////////////
// Public member functions //
/////////////////////////////

void ReplicaExchangeWangLandau::run()
{

  char fileName[51];

  //double currentTime, lastBackUpTime;
  currentTime = lastBackUpTime = MPI_Wtime();
  if (GlobalComm.thisMPIrank == 0)
    printf("Running ReplicaExchangeWangLandau...\n");

  acceptMove = h.checkEnergyInRange(physical_system -> observables[0]);

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
  if (PhysicalSystemComm.thisMPIrank == 0) {
    sprintf(fileName, "energyLatticePos_walker%05d.dat", REWLComm.thisMPIrank);
    physical_system -> writeConfiguration(0, fileName);
    printf("Walker %05d: Start energy = %6.3f\n", REWLComm.thisMPIrank, physical_system -> observables[0]);
  }


//-------------- End initialization --------------//

// WL procedure starts here
  //if (PhysicalSystemComm.thisMPIrank == 0)
  //  MPI_Allreduce(&(h.modFactor), &MaxModFactor, 1, MPI_DOUBLE, MPI_MAX, REWLComm.communicator);
  getMaxModFactor();

  //while (h.modFactor > h.modFactorFinal) {
  while (MaxModFactor > h.modFactorFinal) {
    h.histogramFlat = false;

    while (!(h.histogramFlat)) {
      for (unsigned int MCSteps=0; MCSteps<h.histogramCheckInterval; MCSteps++) {

      //========== One MC move update ========== //

        physical_system -> doMCMove();
        physical_system -> getObservables();

        // check if the energy falls within the energy range
        if ( h.checkEnergyInRange(physical_system -> observables[0]) ) {
          // determine WL acceptance
          if ( exp(h.getDOS(physical_system -> oldObservables[0]) - 
                   h.getDOS(physical_system -> observables[0])) > getRandomNumber2() )
            acceptMove = true;
          else
            acceptMove = false;
        }
        else        
          acceptMove = false;

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

      //========== One MC move update ========== //

      //====== One Replica-exchange update ======//

        if (MCSteps % replicaExchangeInterval == 0) {
          if (replicaExchange()) {
            physical_system -> getObservables();
            h.updateHistogramDOS(physical_system -> observables[0]);
            h.acceptedMoves++;
            physical_system -> acceptMCMove();
          }
        }

      //====== One Replica-exchange update ======//

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
      h.totalMCsteps += h.histogramCheckInterval;
      //h.writeHistogramDOSFile("hist_dos_checkpoint.dat");

      // Check histogram flatness
      h.histogramFlat = h.checkHistogramFlatness();

      // Refresh histogram if needed
      //if (h.numHistogramNotImproved >= h.histogramRefreshInterval)
      //  h.refreshHistogram();

      // Get the maximum ModFactor among all processors
      if (PhysicalSystemComm.thisMPIrank == 0)
        MPI_Allreduce(&(h.modFactor), &MaxModFactor, 1, MPI_DOUBLE, MPI_MAX, REWLComm.communicator);

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

//////////////////////////////
// Private member functions //
//////////////////////////////

bool ReplicaExchangeWangLandau::replicaExchange()
{

  std::cout << "YingWai's check: Inside replicaExchange\n";

  double localDOSRatio             {0.0};
  bool   replicaExchangeAcceptance {false};
  ObservableType energyForExchange = physical_system -> observables[0];
  partnerID = -1;

  // Everyone finds its swap-partner
  assignSwapPartner();
  printf("GlobalID %05d, Walker %05d: partnerID = %05d\n", GlobalComm.thisMPIrank, REWLComm.thisMPIrank, partnerID);

  if ((partnerID != -1) && (PhysicalSystemComm.thisMPIrank == 0)) {
    // Exchange energy with partner
    exchangeEnergy(energyForExchange);

    // If the energy received is within my energy range,
    // calculate DOS ratio from my histogram
    if ( h.checkEnergyInRange(energyForExchange) )
      localDOSRatio = exp( h.getDOS(physical_system -> observables[0]) - h.getDOS(energyForExchange) );
    else localDOSRatio = 0.0;
    
    // Exchange DOS ratios and calculate acceptance probability
    replicaExchangeAcceptance = determineAcceptance(localDOSRatio);

    // Exchange configurations
    if (replicaExchangeAcceptance) {

      physical_system -> observables[0] = energyForExchange;
      exchangeConfiguration(physical_system -> pointerToConfiguration, 1, physical_system -> MPI_ConfigurationType );

      // Up to this point, energy = physical_system -> observables[0] and the configuration are new
      // 1. other observables still need to be calculated, if needed
      // 2. calculate energy from scratch next time

    }
    else {
      //what to do when a replica exchange move is rejected?
    }

    // Performance statistics
  
  }

  //REWLComm.barrier();

  return replicaExchangeAcceptance;

}


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


void ReplicaExchangeWangLandau::exchangeEnergy(ObservableType &energyForSwap)
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


void ReplicaExchangeWangLandau::exchangeConfiguration(void* ptrToConfig, int numElements, MPI_Datatype MPI_config_type)
{
  if (partnerID != -1) {
    REWLComm.swapVector(ptrToConfig, numElements, MPI_config_type, partnerID);
    physical_system -> getObservablesFromScratch = true;
  }

}


void ReplicaExchangeWangLandau::getMaxModFactor()
{

  //YingWai's note: will it be more performant if it is split into two steps?  (Dec 25, 17)
  // * Asynchronized AllReduce within REWLComm
  // * Broadcast to group members within PhySystemComm
  MPI_Allreduce(&(h.modFactor), &MaxModFactor, 1, MPI_DOUBLE, MPI_MAX, GlobalComm.communicator);

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
          if (key == "replicaExchangeInterval") {
            lineStream >> replicaExchangeInterval;
            //std::cout << "REWL: replicaExchangeInterval = " << replicaExchangeInterval << std::endl;
            continue;
          }
  
        }

      }

    }
    inputFile.close();
  }


}

