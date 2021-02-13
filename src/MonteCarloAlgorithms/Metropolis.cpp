#include <cassert>
#include <cmath>
#include <filesystem>
#include <sstream>
#include "Metropolis.hpp"
#include "Utilities/RandomNumberGenerator.hpp"
#include "Utilities/CheckFile.hpp"
#include "Utilities/CompareNumbers.hpp"

// Constructor
Metropolis::Metropolis(PhysicalSystem* ps, const char* inputFile)
{

  if (GlobalComm.thisMPIrank == 0)
    printf("\nSimulation method: Metropolis sampling \n");

  if (std::filesystem::exists(inputFile))
    readMCInputFile(inputFile);
  else {
    std::cout << "Error: No input file for reading Metropolis simulation info. Quiting... \n";
    exit(7);
  }

  physical_system = ps;

  // Allocate space to store observables and other statistics
  if (physical_system->numObservables > 0) {

     averagedObservables         = new ObservableType[physical_system->numObservables];
     averagedObservablesSquared  = new ObservableType[physical_system->numObservables];
     standardDeviations          = new ObservableType[physical_system->numObservables];

    for (unsigned int i=0; i<physical_system->numObservables; i++) {
      averagedObservables[i]         = 0.0;
      averagedObservablesSquared[i]  = 0.0;
      standardDeviations[i]          = 0.0;
    }

  }

  if (std::filesystem::exists("mc.dat")) 
    timeSeriesFile = fopen("mc.dat", "a");
  else 
    timeSeriesFile = fopen("mc.dat", "w"); 

  fprintf(timeSeriesFile, "# Thermalization: (%lu steps) \n", numberOfThermalizationSteps);
  fprintf(timeSeriesFile, "# Temperature %8.5f\n", temperature);
  fprintf(timeSeriesFile, "# MC steps           Observables\n");

  if (simInfo.restartFlag) {
    if (std::filesystem::exists("metropolis_checkpoint.dat"))
      readCheckPointFile("metropolis_checkpoint.dat");
    else {
      std::cout << "\n   WARNING! Restart file 'metropolis_checkpoint.dat' not found. ";
      std::cout << "\n            Performing a fresh run instead of a restarted run. \n\n";
    }
  }

}


//Destructor
Metropolis::~Metropolis()
{

  delete[] averagedObservables;
  delete[] averagedObservablesSquared;
  delete[] standardDeviations;

  fclose(timeSeriesFile);

  if (GlobalComm.thisMPIrank == 0)
    printf("Exiting Metropolis class... \n");

}


void Metropolis::run() 
{

  char fileName[51];

  currentTime = lastBackUpTime = MPI_Wtime();
  if (GlobalComm.thisMPIrank == 0)
    printf("   Running Metropolis Sampling...\n");

  // Thermalization (observables are not accumulated)
  while (thermalizationStepsPerformed < numberOfThermalizationSteps) {
  //for (unsigned long int MCSteps=0; MCSteps<numberOfThermalizationSteps; MCSteps++) {

    for (unsigned long int i=0; i<numberOfMCUpdatesPerStep; i++) {

      physical_system -> doMCMove();
      physical_system -> getObservables();
  
      // Determine acceptance
      if ( exp((physical_system -> oldObservables[0] - physical_system -> observables[0]) / temperature ) > getRandomNumber2() )
        physical_system -> acceptMCMove(); 
      else
        physical_system -> rejectMCMove();

    }

    thermalizationStepsPerformed++;
    physical_system -> getAdditionalObservables();       // can skip this to save time

    // Write observables to file
    writeMCFile(thermalizationStepsPerformed);
    
    // Write restart files at interval
    currentTime = MPI_Wtime();
    if (GlobalComm.thisMPIrank == 0) {
      if (currentTime - lastBackUpTime > checkPointInterval) {
        writeCheckPointFiles(checkPoint);
        lastBackUpTime = currentTime;
      }
    }
 
  }

  fprintf(timeSeriesFile, "# End of thermalization. \n\n");
  fprintf(timeSeriesFile, "# Accumulation: (%lu steps) \n", numberOfMCSteps);
  fprintf(timeSeriesFile, "# Temperature %8.5f\n", temperature);
  fprintf(timeSeriesFile, "# MC steps           Observables\n");

  // Observable accumulation starts here
  while (MCStepsPerformed < numberOfMCSteps) {
  //for (unsigned long int MCSteps=0; MCSteps<numberOfMCSteps; MCSteps++) {

    for (unsigned long int i=0; i<numberOfMCUpdatesPerStep; i++) {
    
      physical_system -> doMCMove();
      physical_system -> getObservables();
  
      // Determine acceptance
      if ( exp((physical_system -> oldObservables[0] - physical_system -> observables[0]) / temperature ) > getRandomNumber2() ) {
        physical_system -> acceptMCMove();
        acceptedMoves++;
      }
      else {
        physical_system -> rejectMCMove();
        rejectedMoves++;
      }

    }
    MCStepsPerformed++;

    physical_system -> getAdditionalObservables();
    accumulateObservables(); 

    // Write observables to file
    writeMCFile(MCStepsPerformed);

    // Write restart files at interval
    currentTime = MPI_Wtime();
    if (GlobalComm.thisMPIrank == 0) {

      if (currentTime - lastBackUpTime > checkPointInterval) {
        writeCheckPointFiles(checkPoint);
        lastBackUpTime = currentTime;
      }

      if (MCStepsPerformed % configurationWriteInterval == 0) {
        sprintf(fileName, "configurations/config%012lu.dat", MCStepsPerformed);
        physical_system -> writeConfiguration(0, fileName);
        sprintf(fileName, "configurations/config%012lu.xyz", MCStepsPerformed);
        physical_system -> writeConfiguration(1, fileName);
      }

    }

  }
  fprintf(timeSeriesFile, "# End of accumulation. \n\n");
  writeCheckPointFiles(checkPoint);

  calculateAveragesAndVariances();
  writeCheckPointFiles(endOfSimulation);
  
  physical_system -> calculateThermodynamics(averagedObservables, averagedObservablesSquared, temperature);

}


// This implementation is similar to the one in Histogram class
void Metropolis::readMCInputFile(const char* fileName)
{

  if (GlobalComm.thisMPIrank == 0) 
    std::cout << "   Metropolis class reading input file: " << fileName << "\n";

  std::ifstream inputFile(fileName);   // TODO: check if a file stream is initialized
  std::string line, key;

  if (inputFile.is_open()) {
    
    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
        
        std::istringstream lineStream(line);
        lineStream >> key;

        if (key.compare(0, 1, "#") != 0) {

          if (key == "numberOfThermalizationSteps") {
            lineStream >> numberOfThermalizationSteps;
            //std::cout << "Metropolis: numberOfThermalizationSteps = " << numberOfThermalizationSteps << "\n";
            continue;
          }
          else if (key == "numberOfMCSteps") {
            lineStream >> numberOfMCSteps;
            //std::cout << "Metropolis: numberOfMCSteps = " << numberOfMCSteps << "\n";
            continue;
          }
          else if (key == "numberOfMCUpdatesPerStep") {
            lineStream >> numberOfMCUpdatesPerStep;
            //std::cout << "Metropolis: numberOfMCUpdatesPerStep = " << numberOfMCUpdatesPerStep << "\n";
            continue;
          }
          else if (key == "temperature") {
            lineStream >> temperature;
            //std::cout << "Metropolis: temperature = " << temperature << "\n";
            continue;
          }
          else if (key == "checkPointInterval") {
            lineStream >> checkPointInterval;
            //std::cout << "Metropolis: checkPointInterval = " << checkPointInterval << " seconds \n";
            continue;
          }
          else if (key == "configurationWriteInterval") {
            lineStream >> configurationWriteInterval;
            //std::cout << "Metropolis: configurationWriteInterval = " << configurationWriteInterval << " seconds \n";
            continue;
          }
          

        }

      }

    }
    inputFile.close();

  }

}


void Metropolis::readCheckPointFile(const char* fileName)
{

  if (GlobalComm.thisMPIrank == 0) 
    std::cout << "   Metropolis class reading checkpoint file: " << fileName << "\n";

  std::ifstream inputFile(fileName);   // TODO: check if a file stream is initialized
  std::string line, key;

  if (inputFile.is_open()) {
    
    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
        
        std::istringstream lineStream(line);
        lineStream >> key;

        if (key.compare(0, 1, "#") != 0) {

          if (key == "temperature") {
            lineStream >> restartTemperature;
            //std::cout << "Metropolis: restartTemperature = " << restartTemperature << "\n";
            continue;
          }
          else if (key == "thermalizationStepsPerformed") {
            lineStream >> thermalizationStepsPerformed;
            //std::cout << "Metropolis: thermalizationStepsPerformed = " << thermalizationStepsPerformed << "\n";
            continue;
          }
          else if (key == "MCStepsPerformed") {
            lineStream >> MCStepsPerformed;
            //std::cout << "Metropolis: MCStepsPerformed = " << MCStepsPerformed << "\n";
            continue;
          }
          else if (key == "acceptedMoves") {
            lineStream >> acceptedMoves;
            //std::cout << "Metropolis: acceptedMoves = " << acceptedMoves << "\n";
            continue;
          }
          else if (key == "rejectedMoves") {
            lineStream >> rejectedMoves;
            //std::cout << "Metropolis: rejectedMoves = " << rejectedMoves << "\n";
            continue;
          }
          else if (key == "averagedObservables") {
            unsigned int counter = 0;
            while (lineStream && counter < physical_system->numObservables) {
              lineStream >> averagedObservables[counter];
              //std::cout << "Metropolis: averageObservables[" << counter << "] = " << averagedObservables[counter] << "\n";
              counter++;
            }
            continue;
          }
          else if (key == "standardDeviations") {
            unsigned int counter = 0;
            while (lineStream && counter < physical_system->numObservables) {
              lineStream >> standardDeviations[counter];
              //std::cout << "Metropolis: standardDeviations[" << counter << "] = " << standardDeviations[counter] << "\n";
              counter++;
            }
            continue;
          }

        }

      }

    }
    inputFile.close();

  }
  
  // Check consistency: temperature
  if (!sameMagnitude(restartTemperature, temperature)) {
    printf("\n   CAUTION! Temperature of previous run different from this run:");
    printf("\n            - Temperature in main input file: %8.5f \n", temperature);
    printf("\n            - Temperature in checkpoint file: %8.5f \n\n", restartTemperature);
    printf("\n            No further work will be performed. Quitting OWL...\n\n");
    exit(7);
  }
  
  // Check consistency: thermalizationSteps
  if (thermalizationStepsPerformed >= numberOfThermalizationSteps) {
    std::cout << "\n   CAUTION! Thermalization steps performed from previous run >= numberOfThermalizationSteps in this run.";
    if (std::filesystem::exists("configurations/config_checkpoint.dat"))
      std::cout << "\n            - Configuration checkpoint file is found, thermalization will be skipped.\n\n";
    else {          // perform thermalization
      thermalizationStepsPerformed = 0;
      std::cout << "\n            - Configuration checkpoint file is not found, thermalization will be performed.\n\n";
    }    
  }

  // Check consistency: MCSteps
  if (MCStepsPerformed >= numberOfMCSteps) {
    std::cout << "\n   CAUTION! Number of MC steps performed from previous run >= numberOfMCSteps required.";
    std::cout << "\n            No further work will be performed. Quitting OWL...\n\n";
    exit(7);
  }

  // Check consistency: acceptedMoves and rejectedMoves
  assert (acceptedMoves + rejectedMoves == MCStepsPerformed * numberOfMCUpdatesPerStep);

  // Restore averagedObservables and averagedObservablesSquared for accumulation
  for (unsigned int i=0; i<physical_system->numObservables; i++) {
    averagedObservablesSquared[i] = (standardDeviations[i] * standardDeviations[i] + averagedObservables[i] * averagedObservables[i]) * double(MCStepsPerformed);
    averagedObservables[i] *= double(MCStepsPerformed);
  }
    
}


void Metropolis::accumulateObservables()
{

  for (unsigned int i=0; i<physical_system->numObservables; i++) {
    averagedObservables[i] += physical_system -> observables[i];
    averagedObservablesSquared[i]  += physical_system -> observables[i] * physical_system -> observables[i];
  }

}


void Metropolis::calculateAveragesAndVariances()
{

  for (unsigned int i=0; i<physical_system->numObservables; i++) {
    averagedObservables[i] /= double(numberOfMCSteps);
    averagedObservablesSquared[i] /= double(numberOfMCSteps);
    standardDeviations[i] = sqrt( averagedObservablesSquared[i] - averagedObservables[i] * averagedObservables[i] );
    //standardDeviations[i] = sqrt((averagedObservablesSquared[i] - averagedObservables[i] * averagedObservables[i]) / double(numberOfMCSteps - 1));
  }

}


void Metropolis::writeMCFile(unsigned long int MCSteps)
{

  fprintf(timeSeriesFile, "%15lu ", MCSteps);
  
  for (unsigned int i=0; i<physical_system->numObservables; i++)
    fprintf(timeSeriesFile, "%15.6f ", physical_system -> observables[i]);

  fprintf(timeSeriesFile, "\n");

}


void Metropolis::writeStatistics(OutputMode output_mode, const char* filename) 
{
  
  FILE* checkPointFile;
  if (filename != NULL)
    checkPointFile = fopen(filename, "w");
  else checkPointFile = stdout;

  switch(output_mode) {

    case endOfSimulation :

      fprintf(checkPointFile, "\n");
      fprintf(checkPointFile, "             Statistics of Metropolis sampling \n");
      fprintf(checkPointFile, "   ----------------------------------------------------- \n");
      fprintf(checkPointFile, "   Simulation temperature         : %8.5f \n", temperature);
      fprintf(checkPointFile, "   Number of thermalization steps :  %lu \n",  numberOfThermalizationSteps);
      fprintf(checkPointFile, "   Total number of MC steps       :  %lu \n",  numberOfMCSteps);
      fprintf(checkPointFile, "   Number of MC updates per step  :  %lu \n",  numberOfMCUpdatesPerStep);
      fprintf(checkPointFile, "   Number of accepted MC updates  :  %lu (%5.2f %%) \n", 
              acceptedMoves, double(acceptedMoves) / double(numberOfMCSteps * numberOfMCUpdatesPerStep) * 100.0);
      fprintf(checkPointFile, "   Number of rejected MC updates  :  %lu (%5.2f %%) \n", 
              rejectedMoves, double(rejectedMoves) / double(numberOfMCSteps * numberOfMCUpdatesPerStep) * 100.0);
      
      fprintf(checkPointFile, "\n");
    
      fprintf(checkPointFile, "                        Observable                          Mean          Std. deviation \n");
      fprintf(checkPointFile, "   --------------------------------------------------------------------------------------- \n");
      for (unsigned int i=0; i<physical_system -> numObservables; i++)
        fprintf(checkPointFile, "   %45s :     %12.5f      %12.5f \n", physical_system -> observableName[i].c_str(), averagedObservables[i], standardDeviations[i]);    
      fprintf(checkPointFile, "\n"); 

      break;

    case checkPoint :

      fprintf(checkPointFile, "temperature                   %8.5f\n", temperature);
      fprintf(checkPointFile, "thermalizationStepsPerformed   %lu\n",  thermalizationStepsPerformed);
      fprintf(checkPointFile, "MCStepsPerformed               %lu\n",  MCStepsPerformed);
      fprintf(checkPointFile, "acceptedMoves                  %lu\n",  acceptedMoves);
      fprintf(checkPointFile, "rejectedMoves                  %lu\n",  rejectedMoves);

      fprintf(checkPointFile, "averagedObservables   ");
      for (unsigned int i=0; i<physical_system -> numObservables; i++)
        fprintf(checkPointFile, "%12.5f      ", averagedObservables[i] / double(MCStepsPerformed));
      fprintf(checkPointFile, "\n");
      
      fprintf(checkPointFile, "standardDeviations   ");
      for (unsigned int i=0; i<physical_system -> numObservables; i++) {
        double temp_ave  = averagedObservables[i] / double(MCStepsPerformed);
        double temp_ave2 = averagedObservablesSquared[i] / double(MCStepsPerformed);
        standardDeviations[i] = sqrt(temp_ave2 - temp_ave*temp_ave);
        fprintf(checkPointFile, "%12.5f      ", standardDeviations[i]);
      }
      fprintf(checkPointFile, "\n");
      
      break;

    default :
      break;
      
  }

  if (filename != NULL) fclose(checkPointFile);

}


void Metropolis::writeCheckPointFiles(OutputMode output_mode)
{

  char fileName[51];

  switch (output_mode) {

    case endOfIteration :
      break;

    case endOfSimulation :
      sprintf(fileName, "configurations/config_final.dat");
      physical_system -> writeConfiguration(0, fileName);
      writeStatistics(endOfSimulation);
      writeStatistics(endOfSimulation, "metropolis_final.dat");
      break;

    case checkPoint :
      physical_system -> writeConfiguration(0, "configurations/config_checkpoint.dat");
      writeStatistics(checkPoint, "metropolis_checkpoint.dat"); 
      break;

    default :
      break;

  };


}
