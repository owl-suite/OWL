//#include <iostream>
#include "Communications.hpp"
#include "RandomNumberGenerator.hpp"
#include "InputOutput.hpp"
#include "OWLMain.hpp"


SimulationInfo simInfo;

int main (int argc, char *argv[]) {

  PhysicalSystem*       physical_system;
  MonteCarloAlgorithm*  MC;
  MPICommunicator       physicalSystemComm;
  MPICommunicator       mcAlgorithmComm;


  readCommandLineArguments(argc, argv);

  // Read from input file; store the simulation settings in simInfo
  readMainInputFile(argv[1], simInfo);

  // Initializations
  initializeMPICommunication(physicalSystemComm, mcAlgorithmComm);
  initializeRandomNumberGenerator(physicalSystemComm, simInfo.rngSeed);

  // Determine and setup the PhysicalSystem and MonteCarlo classes
  setSimulation(physical_system, MC, physicalSystemComm, mcAlgorithmComm);

  // Main calculations
  MC -> run(physical_system); 

  // Finalize
  finalizeMPICommunication();

  delete physical_system;
  delete MC;

  std::cout << "\nSimulation finished :)\n";
  return 0;
}
