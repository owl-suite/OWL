#include <cstdio>
#include <cmath>
#include <filesystem>
#include "GlobalUpdateWangLandauSampling.hpp"
#include "Utilities/RandomNumberGenerator.hpp"
//#include "Communications.hpp"


// Constructor
GlobalUpdateWangLandauSampling::GlobalUpdateWangLandauSampling(PhysicalSystem* ps) : h(simInfo.restartFlag, simInfo.MCInputFile, simInfo.HistogramCheckpointFile)
{

  if (GlobalComm.thisMPIrank == 0)
    printf("\nStarting Global-update Wang-Landau sampling...\n");

  physical_system = ps;

}


//Destructor
GlobalUpdateWangLandauSampling::~GlobalUpdateWangLandauSampling()
{

  if (GlobalComm.thisMPIrank == 0)
    printf("Exiting GlobalUpdateWangLandauSampling class... \n");

}


void GlobalUpdateWangLandauSampling::run()
{

  char fileName[51];

  currentTime = lastBackUpTime = MPI_Wtime();
  if (GlobalComm.thisMPIrank == 0)
    printf("   Finding initial configuration within energy range... ");

  // Find the first energy that falls within the WL energy range    
  while (!acceptMove) {
    physical_system -> doMCMove();
    physical_system -> getObservables();
    physical_system -> acceptMCMove();    // always accept the move to push the state forward
    acceptMove = h.checkEnergyInRange(physical_system -> observables[0]);
  }

  // Always accept the first energy if it is within range
  h.globalUpdateHistogramDOS(physical_system -> observables[0]);

  // Write out the energy
  if (GlobalComm.thisMPIrank == 0 && std::filesystem::exists("configurations")) {
    physical_system -> writeConfiguration(0, "configurations/config_initial.dat");
    printf("done.\n");
  }

//-------------- End initialization --------------//

// WL procedure starts here
  while (h.modFactor > h.modFactorFinal) {
    h.histogramFlat = false;

    if (GlobalComm.thisMPIrank == 0) {
      printf("   Running iteration %2d   (f = %12.8e) ... ", h.iterations, h.modFactor);
      fflush(stdout);
    }

    h.numberOfUpdatesPerIteration = 0;
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
           h.globalUpdateHistogramDOS(physical_system -> observables[0]);
           h.acceptedMoves++;        
 
           physical_system -> acceptMCMove();
        }
        else {
           physical_system -> rejectMCMove();

           // Update histogram and DOS with oldEnergy
           h.globalUpdateHistogramDOS(physical_system -> oldObservables[0]);
           h.rejectedMoves++;
        }
        h.totalMCsteps++;
     
        // Write restart files at interval
        currentTime = MPI_Wtime();
        if (GlobalComm.thisMPIrank == 0) {
          if (currentTime - lastBackUpTime > checkPointInterval) {
            h.writeHistogramDOSFile("hist_dos_checkpoint.dat");
            physical_system -> writeConfiguration(1, "configurations/config_checkpoint.dat");
            lastBackUpTime = currentTime;
          }
        }
      }
      //h.writeHistogramDOSFile("hist_dos_checkpoint.dat");
      h.numberOfUpdatesPerIteration += h.histogramCheckInterval;

      // Check histogram flatness
      h.histogramFlat = h.checkHistogramFlatness();
      
      //if (h.numHistogramNotImproved >= h.histogramRefreshInterval)
      //  h.refreshHistogram();
    }

    //bool KB  = h.checkKullbackLeiblerDivergence();

    if (GlobalComm.thisMPIrank == 0) {
      printf("done.\n");

      // Also write restart files here 
      sprintf(fileName, "hist_dos_iteration%02d.dat", h.iterations);
      h.writeHistogramDOSFile(fileName);
      physical_system -> writeConfiguration(1, "configurations/config_checkpoint.dat");
    }

    // Go to next iteration
    h.modFactor /= h.modFactorReducer;
    h.resetHistogram();
    h.iterations++;
  }

  // Write out data at the end of the simulation
  h.writeNormDOSFile("dos.dat");
  h.writeHistogramDOSFile("hist_dos_checkpoint.dat");
  h.writeHistogramDOSFile("hist_dos_final.dat");

}

