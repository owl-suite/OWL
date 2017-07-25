#include <iostream>
#include "Communications.hpp"
#include "RandomNumberGenerator.hpp"
#include "InputOutput.hpp"
#include "OWLMain.hpp"


int main (int argc, char *argv[]) {

  SimulationInfo        simInfo;
  PhysicalSystem*       physical_system;
  MonteCarloAlgorithm*  MC;
  MPICommunicator       physicalSystemComm;
  MPICommunicator       mcAlgorithmComm;


  readCommandLineArguments(argc, argv);

  // Read from input file
  readMainInputFile(argv[1], simInfo);

  // Initializations
  initializeMPICommunication(simInfo, physicalSystemComm, mcAlgorithmComm);
  initializeRandomNumberGenerator(physicalSystemComm, simInfo.rngSeed);

  // Determine and setup the PhysicalSystem and MonteCarlo classes
  setSimulation(simInfo, physical_system, MC, physicalSystemComm, mcAlgorithmComm);

  // Main calculations
  MC -> run(physical_system); 

  // Finalize
  finalizeMPICommunication(simInfo);

  delete physical_system;
  delete MC;

  std::cout << "\nSimulation finished :)\n";
  return 0;
}
