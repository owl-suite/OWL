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
    printf("Simulation method: Metropolis sampling\n");

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

  if (GlobalComm.thisMPIrank == 0)
    printf("Exiting Metropolis class... \n");

}

void Metropolis::run()
{

  currentTime = lastBackUpTime = MPI_Wtime();
  if (GlobalComm.thisMPIrank == 0)
    printf("Running Metropolis Sampling...\n");

  // Thermalization (observables are not accumulated)
  for (unsigned long int MCSteps=0; MCSteps<numberOfThermalizationSteps; MCSteps++) {

    for (unsigned long int i=0; i<numberOfMCUpdatesPerStep; i++) {

      physical_system -> doMCMove();
      physical_system -> getObservables();
  
      // Determine acceptance
      if ( exp((physical_system -> oldObservables[0] - physical_system -> observables[0]) * temperature ) > getRandomNumber2() )
        physical_system -> acceptMCMove();
      else
        physical_system -> rejectMCMove();
  
    }
 
  }

  // Observable accumulation starts here
  for (unsigned long int MCSteps=0; MCSteps<numberOfMCSteps; MCSteps++) {

    for (unsigned long int i=0; i<numberOfMCUpdatesPerStep; i++) {
    
      physical_system -> doMCMove();
      physical_system -> getObservables();
  
      // Determine acceptance
      if ( exp((physical_system -> oldObservables[0] - physical_system -> observables[0]) * temperature ) > getRandomNumber2() ) {
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

  }

  calculateAveragesAndVariances();

  writeResultsFile();


}

// This implementation is similar to the one in Histogram class
void Metropolis::readMCInputFile(const char* fileName)
{

  if (GlobalComm.thisMPIrank == 0) 
    std::cout << "Metropolis class reading input file: " << fileName << std::endl;

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
          if (key == "numberOfMCSteps") {
            lineStream >> numberOfMCSteps;
            //std::cout << "Metropolis: numberOfMCSteps = " << numberOfMCSteps << std::endl;
            continue;
          }
          if (key == "numberOfMCUpdatesPerStep") {
            lineStream >> numberOfMCUpdatesPerStep;
            //std::cout << "Metropolis: numberOfMCUpdatesPerStep = " << numberOfMCUpdatesPerStep << std::endl;
            continue;
          }
          if (key == "temperature") {
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

  //TODO: print observables to stdout/file

  for (unsigned int i=0; i<physical_system->numObservables; i++) {
    averagedObservables[i] += physical_system -> observables[i];    
    variances[i] += (physical_system -> observables[i]) * (physical_system -> observables[i]);
  }

}

// calculate average and variance of observables
void Metropolis::calculateAveragesAndVariances()
{

  for (unsigned int i=0; i<physical_system->numObservables; i++) {
    averagedObservables[i] /= (double)numberOfMCSteps;
    variances[i] -= averagedObservables[i] * averagedObservables[i];
  }

}


void Metropolis::writeResultsFile(const char* filename)
{

  FILE* f;
  if (filename != NULL) f = fopen(filename, "w");
  else f = stdout;

  fprintf(f, "\n");
  fprintf(f, "Statistics of Metropolis sampling \n");
  fprintf(f, "--------------------------------- \n");
  fprintf(f, "Simulation temperature:         %8.5f \n", temperature);
  fprintf(f, "Total number of MC steps:       %lu \n", numberOfMCSteps);
  fprintf(f, "Number of thermalization steps: %lu \n", numberOfThermalizationSteps);
  fprintf(f, "Number of MC sweeps per steps:  %lu \n", numberOfMCUpdatesPerStep);
  fprintf(f, "Number of accepted MC moves:    %lu (%5.2f %%) \n", 
          acceptedMoves, (double)acceptedMoves / (double)numberOfMCSteps * 100.0);
  fprintf(f, "Number of rejected MC moves:    %lu (%5.2f %%) \n", 
          rejectedMoves, (double)rejectedMoves / (double)numberOfMCSteps * 100.0);
  
  fprintf(f, "Observable         Variance \n");
  fprintf(f, "--------------------------- \n");
  for (unsigned int i=0; i<physical_system -> numObservables; i++)
    fprintf(f, "%10.5f         %10.5f \n", averagedObservables[i], variances[i]);

}