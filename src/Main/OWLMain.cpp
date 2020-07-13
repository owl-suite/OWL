//#include <iostream>
#include "Utilities/RandomNumberGenerator.hpp"
#include "Communications.hpp"
#include "InputOutput.hpp"
#include "OWLMain.hpp"

//SimulationInfo simInfo;

int main (int argc, char *argv[]) {

  PhysicalSystem*       physical_system;
  MonteCarloAlgorithm*  MC;
  MPICommunicator       physicalSystemComm;
  MPICommunicator       mcAlgorithmComm;


  readCommandLineArguments(argc, argv);

  // Read from input file; store the simulation settings in simInfo
  readMainInputFile(argv[1]);

  writeSimulationInfo();

  // Initializations
  initializeMPICommunication( physicalSystemComm, mcAlgorithmComm );
  initializeRandomNumberGenerator( physicalSystemComm, simInfo.rngSeed );

  // Determine and setup the PhysicalSystem and MonteCarlo classes
  setSimulation( physical_system, MC, physicalSystemComm, mcAlgorithmComm );

  // Main calculations
  MC -> run(); 

  // Finalize simulation and MPI
  finalizeSimulation( physical_system, MC );
  finalizeMPICommunication();

  if (GlobalComm.thisMPIrank == 0)
    std::cout << "\nSimulation finished :)\n\n";

  return 0;
}
