#include <iostream>
#include "Communications.hpp"
#include "RandomNumberGenerator.hpp"
#include "InputOutput.hpp"
#include "OWLMain.hpp"


int main (int argc, char *argv[]) {

  PhysicalSystem*      physical_system;
  MonteCarloAlgorithm* MC;
  SimulationInfo       simInfo;

  readCommandLineArguments(argc, argv);

  // Read from input file
  readMainInputFile(argv[1], simInfo);

  // Initializations
  initializeMPICommunication(simInfo);
  initializeRandomNumberGenerator(simInfo.rngSeed);

  // Determine and setup the PhysicalSystem and MonteCarlo classes
  setSimulation(physical_system, MC, simInfo);

  // Main calculations
  MC -> run(physical_system); 

  // Finalize
  finalizeMPICommunication(simInfo);

  delete physical_system;
  delete MC;

  std::cout << "\nSimulation finished :)\n";
  return 0;
}
