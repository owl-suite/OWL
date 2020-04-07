// MUCA implementation following the recipe suggested in this paper:
// Ref: J. Gross, J. Zierenberg, M. Weigel, and W. Janke. Comp. Phys. Comm. 224, 387–395 (2018). 

#include <cstdio>
#include <cmath>
#include "MulticanonicalSampling.hpp"
#include "Utilities/RandomNumberGenerator.hpp"


// Constructor
MulticanonicalSampling::MulticanonicalSampling(PhysicalSystem* ps) : h(simInfo.restartFlag, simInfo.MCInputFile, simInfo.HistogramCheckpointFile)
{

  printf("Simulation method: Multicanonical (MUDA) sampling\n");
  physical_system = ps;

}


// Destructor
MulticanonicalSampling::~MulticanonicalSampling()
{

  printf("Exiting MulticanonicalSampling class... \n");

}


void MulticanonicalSampling::run()
{

  currentTime = lastBackUpTime = MPI_Wtime();
  if (GlobalComm.thisMPIrank == 0)
    printf("Running Multicanonical Sampling...\n");

  char fileName[51];
 
  //-------------- Initialization starts --------------//

  // Find the first energy that falls within the energy range
  while (!acceptMove) {
    physical_system -> doMCMove();
    physical_system -> getObservables();
    physical_system -> acceptMCMove();    // always accept the move to push the state forward
    acceptMove = h.checkEnergyInRange(physical_system -> observables[0]);
  }

  // Write out the initial configuration
  if (GlobalComm.thisMPIrank == 0)
    physical_system -> writeConfiguration(0, "initial_configuration.dat");

//---------------- Initialization ends ----------------//

  // MUCA procedure starts here
  while (!(h.histogramFlat)) {

    // Thermalization (these steps do not update the histogram)
    for (unsigned int MCSteps=0; MCSteps<h.numberOfThermalizationSteps; MCSteps++) {

      physical_system -> doMCMove();
      physical_system -> getObservables();

      // check if the energy falls within the energy range
      if ( !h.checkEnergyInRange(physical_system -> observables[0]) ) {
        acceptMove = false;
        physical_system -> rejectMCMove();
      }
      else {
        // determine acceptance
        if ( exp(h.getDOS(physical_system -> oldObservables[0]) - 
                 h.getDOS(physical_system -> observables[0])) > getRandomNumber2() ) {
          acceptMove = true;
          physical_system -> acceptMCMove();
                 }
        else {
          acceptMove = false;
          physical_system -> rejectMCMove();
        }
      }

    }
    h.totalMCsteps += h.numberOfThermalizationSteps;

    // MUCA statistics starts here
    for (unsigned int MCSteps=0; MCSteps<h.numberOfUpdatesPerIteration; MCSteps++) {

      physical_system -> doMCMove();
      physical_system -> getObservables();

      // check if the energy falls within the energy range
      if ( !h.checkEnergyInRange(physical_system -> observables[0]) )
        acceptMove = false;
      else {
        // determine acceptance
        if ( exp(h.getDOS(physical_system -> oldObservables[0]) - 
                 h.getDOS(physical_system -> observables[0])) > getRandomNumber2() )
          acceptMove = true;
        else
          acceptMove = false;
      }

      if (acceptMove) {
         // Update histogram with trial state
         h.updateHistogram(physical_system -> observables[0]);
         h.acceptedMoves++;        

         physical_system -> acceptMCMove();
      }
      else {
         physical_system -> rejectMCMove();

         // Update histogram with old state
         h.updateHistogram(physical_system -> oldObservables[0]);
         h.rejectedMoves++;
      }
   
      // Write restart files at interval
      currentTime = MPI_Wtime();
      if (currentTime - lastBackUpTime > 300) {
        if (GlobalComm.thisMPIrank == 0) {
          h.writeHistogramDOSFile("hist_dos_checkpoint.dat");
          physical_system -> writeConfiguration(1, "OWL_restart_configuration");
          lastBackUpTime = currentTime;
        }
      }
    }
    h.totalMCsteps += h.numberOfUpdatesPerIteration;

    // Update DOS with the histogram
    h.updateDOSwithHistogram();
    // check flatness of histogram using Kullback-Leibler divergence
    h.histogramFlat = h.checkKullbackLeiblerDivergence();
      
    if (GlobalComm.thisMPIrank == 0) {
      printf("Iteration %d finished \n", h.iterations);
    
      // Also write restart file here 
      sprintf(fileName, "hist_dos_iteration%02d.dat", h.iterations);
      h.writeHistogramDOSFile(fileName);
      physical_system -> writeConfiguration(1, "OWL_restart_input");
    }

    // Go to next iteration
    h.resetHistogram();
    h.iterations++;
    h.numberOfUpdatesPerIteration = static_cast<unsigned int>(ceil(static_cast<double>(h.numberOfUpdatesPerIteration) * h.numberOfUpdatesMultiplier));
    if (GlobalComm.thisMPIrank == 0)
      printf("Number of updates in the next iteration = %d\n \n", h.numberOfUpdatesPerIteration);

  }

  // Write out data at the end of the simulation
  if (GlobalComm.thisMPIrank == 0) {
    h.writeNormDOSFile("dos.dat");
    h.writeHistogramDOSFile("hist_dos_final.dat");
    printf("Number of total MC steps (including thermalization) = %lu\n", h.totalMCsteps);
  }

}

