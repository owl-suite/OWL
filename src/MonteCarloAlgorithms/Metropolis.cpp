#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include "Metropolis.hpp"
#include "Utilities/RandomNumberGenerator.hpp"
#include "Utilities/CheckFile.hpp"

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

  // Allocate space to store observables and their squares for calculating variances
  if (physical_system -> numObservables > 0) {
     averagedObservables = new ObservableType[physical_system->numObservables];
     variances = new ObservableType[physical_system->numObservables];
  }

  for (unsigned int i=0; i<physical_system->numObservables; i++) {
    averagedObservables[i] = 0.0;
    variances[i] = 0.0;
  }

  if (!std::filesystem::exists("configurations"))
    std::filesystem::create_directory("configurations");

  if (std::filesystem::exists("mc.dat"))
    MCOutputFile = fopen("mc.dat", "a");
  else {
    MCOutputFile = fopen("mc.dat", "w");
    fprintf(MCOutputFile, "# MC steps           Observables\n");
  }

}


//Destructor
Metropolis::~Metropolis()
{

  delete[] averagedObservables;
  delete[] variances;

  fclose(MCOutputFile);

  if (GlobalComm.thisMPIrank == 0)
    printf("Exiting Metropolis class... \n");

}


void Metropolis::run()
{

  currentTime = lastBackUpTime = MPI_Wtime();
  if (GlobalComm.thisMPIrank == 0)
    printf("   Running Metropolis Sampling...\n");

  char fileName[51];

  // Thermalization (observables are not accumulated)
  for (unsigned long int MCSteps=0; MCSteps<numberOfThermalizationSteps; MCSteps++) {

    for (unsigned long int i=0; i<numberOfMCUpdatesPerStep; i++) {

      physical_system -> doMCMove();
      physical_system -> getObservables();
  
      // Determine acceptance
      if ( exp((physical_system -> oldObservables[0] - physical_system -> observables[0]) / temperature ) > getRandomNumber2() )
        physical_system -> acceptMCMove(); 
      else
        physical_system -> rejectMCMove();

    }

    // Write observables to file
    writeMCFile(MCSteps);
    
    // Write restart files at interval
    currentTime = MPI_Wtime();
    if (GlobalComm.thisMPIrank == 0) {
      if (currentTime - lastBackUpTime > checkPointInterval) {
        physical_system -> writeConfiguration(1, "config_checkpoint.dat");
        lastBackUpTime = currentTime;
      }
    }
 
  }

  // Observable accumulation starts here
  for (unsigned long int MCSteps=0; MCSteps<numberOfMCSteps; MCSteps++) {

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
    
    // Accumulate observables
    accumulateObservables(); 

    // Write observables to file
    writeMCFile(MCSteps);

    // Write restart files at interval
    currentTime = MPI_Wtime();
    if (GlobalComm.thisMPIrank == 0) {

      if (currentTime - lastBackUpTime > checkPointInterval) {
        physical_system -> writeConfiguration(1, "config_checkpoint.dat");
        lastBackUpTime = currentTime;
      }

      if (MCSteps % configurationWriteInterval == 0) {
        sprintf(fileName, "configurations/config%012lu.dat", MCSteps);
        physical_system -> writeConfiguration(1, fileName);
      } 

    }

  }

  calculateAveragesAndVariances();

  writeResultsFile();

}


// This implementation is similar to the one in Histogram class
void Metropolis::readMCInputFile(const char* fileName)
{

  if (GlobalComm.thisMPIrank == 0) 
    std::cout << "   Metropolis class reading input file: " << fileName << std::endl;

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
            //std::cout << "Metropolis: numberOfThermalizationSteps = " << numberOfThermalizationSteps << std::endl;
            continue;
          }
          else if (key == "numberOfMCSteps") {
            lineStream >> numberOfMCSteps;
            //std::cout << "Metropolis: numberOfMCSteps = " << numberOfMCSteps << std::endl;
            continue;
          }
          else if (key == "numberOfMCUpdatesPerStep") {
            lineStream >> numberOfMCUpdatesPerStep;
            //std::cout << "Metropolis: numberOfMCUpdatesPerStep = " << numberOfMCUpdatesPerStep << std::endl;
            continue;
          }
          else if (key == "temperature") {
            lineStream >> temperature;
            //std::cout << "Metropolis: temperature = " << temperature << std::endl;
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


void Metropolis::accumulateObservables()
{

  for (unsigned int i=0; i<physical_system->numObservables; i++) {
    averagedObservables[i] += physical_system -> observables[i];    
    variances[i] += (physical_system -> observables[i]) * (physical_system -> observables[i]);
  }

}


// calculate average and variance of observables
void Metropolis::calculateAveragesAndVariances()
{

  for (unsigned int i=0; i<physical_system->numObservables; i++) {
    averagedObservables[i] /= double(numberOfMCSteps);
    variances[i] /= double(numberOfMCSteps);
    variances[i] = sqrt( variances[i] - averagedObservables[i] * averagedObservables[i] );
    //variances[i] = sqrt((variances[i] - averagedObservables[i] * averagedObservables[i]) / double(numberOfMCSteps - 1));
  }

}


void Metropolis::writeMCFile(unsigned long int MCSteps)
{

  fprintf(MCOutputFile, "%15lu ", MCSteps);
  
  for (unsigned int i=0; i<physical_system->numObservables; i++)
    fprintf(MCOutputFile, "%15.6f ", physical_system -> observables[i]);

  fprintf(MCOutputFile, "\n");

}


void Metropolis::writeResultsFile(const char* filename) 
{

  FILE* resultsFile;
  if (filename != NULL) resultsFile = fopen(filename, "w");
  else resultsFile = stdout;

  fprintf(resultsFile, "\n");
  fprintf(resultsFile, "   Statistics of Metropolis sampling \n");
  fprintf(resultsFile, "   --------------------------------- \n");
  fprintf(resultsFile, "   Simulation temperature         : %8.5f \n", temperature);
  fprintf(resultsFile, "   Total number of MC steps       :  %lu \n", numberOfMCSteps);
  fprintf(resultsFile, "   Number of thermalization steps :  %lu \n", numberOfThermalizationSteps);
  fprintf(resultsFile, "   Number of MC updates per step  :  %lu \n", numberOfMCUpdatesPerStep);
  fprintf(resultsFile, "   Number of accepted MC updates  :  %lu (%5.2f %%) \n", 
          acceptedMoves, double(acceptedMoves) / double(numberOfMCSteps * numberOfMCUpdatesPerStep) * 100.0);
  fprintf(resultsFile, "   Number of rejected MC updates  :  %lu (%5.2f %%) \n", 
          rejectedMoves, double(rejectedMoves) / double(numberOfMCSteps * numberOfMCUpdatesPerStep) * 100.0);
  
  fprintf(resultsFile, "\n");

  fprintf(resultsFile, "                        Observable                          Mean          Std. deviation \n");
  fprintf(resultsFile, "   -------------------------------------------------------------------------------------- \n");
  for (unsigned int i=0; i<physical_system -> numObservables; i++)
    fprintf(resultsFile, "   %45s :     %12.5f      %12.5f \n", physical_system -> observableName[i].c_str(), averagedObservables[i], variances[i]);

  fprintf(resultsFile, "\n");

  if (filename != NULL) fclose(resultsFile);

}