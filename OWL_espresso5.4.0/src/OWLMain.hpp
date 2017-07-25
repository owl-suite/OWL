#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include "Globals.hpp"
#include "MCAlgorithms.hpp"
#include "WangLandauSampling.hpp"
#include "MulticanonicalSampling.hpp"
#include "HistogramFreeMUCA.hpp"
#include "Heisenberg2D.hpp"
#include "Ising2D.hpp"
#include "QuantumEspressoSystem.hpp"


void readCommandLineArguments(int argc, char* argv[]) {

  // Read in restart information
  if (argc != 2) {
    //std::cout << "Usage: " << argv[0] << " [RestartFlag=0/1] \n";
    std::cout << "Usage: " << argv[0] << " [Input File Name] \n";
    exit(7);
  }
/*
  else {
    if (sscanf(argv[1], "%d", &restartFlag) != 1) {
      std::cout << "Error: cannot read restart flag! \n";
      exit(1);
    }
  }

  switch (restartFlag) {
    case 0 : 
      // fresh start
      std::cout << "YingWai's check: Fresh MC simulation\n";
      break;
    case 1 : 
      // restart from file
      std::cout << "YingWai's check: Restarted run\n";
      break;
    default :
      std::cout << "Restart flag != 0 or 1. Assume a fresh start.\n";
      restartFlag = 0;
    }
*/
}

void setSimulation(SimulationInfo       simInfo,
                   PhysicalSystem*      &physical_system,
                   MonteCarloAlgorithm* &MC,
                   MPICommunicator      physicalSystemComm,
                   MPICommunicator      mcAlgorithmComm)
{
  
  // Determine MC algorithm
  // 1. Metropolis sampling
  // 2. Wang-Landau sampling
  // 3. Multicanonical sampling (MUCA)
  // 4. Parallel tempering
  // 5. Replica-Exchange Wang-Landau sampling (REWL)
  // 6. Histogram-free Multicanonical sampling (discrete energy version)
  switch (simInfo.algorithm) {
    case 1 :
      MC = new Metropolis( simInfo.restartFlag, simInfo.MCInputFile );
      break;

    case 2 :
      MC = new WangLandauSampling( simInfo.restartFlag, simInfo.MCInputFile );
      break;

    case 3 :
      MC = new MulticanonicalSampling( simInfo.restartFlag, simInfo.MCInputFile );
      break;

    case 4 :
      std::cout << "Parallel tempering not yet implemented\n";
      exit(10);
      break;

    case 5 :
      std::cout << "Replica-Exchange Wang-Landau (REWL) not yet implemented\n";
      exit(10);
      break;
      
    case 6 :
      MC = new DiscreteHistogramFreeMUCA( simInfo.restartFlag, simInfo.MCInputFile );
      break;
      
    default :
      std::cout << "Monte Carlo algorithm not specified. Use default: Wang-Landau sampling.\n";
  }

  // Physical System #####
  //  1: QuantumExpresso
  //  2: LSMS  
  //  3: Heisenberg 2D
  //  4: Ising 2D
  //  5: ...
  switch (simInfo.system) {
    case 1 :
      physical_system = new QuantumEspressoSystem(simInfo, physicalSystemComm);
      break;

    case 2 :
      std::cout << "LSMS not yet available.\n";
      break;

    case 3 :
      physical_system = new Heisenberg2D(simInfo);
      break;

    case 4 :
      physical_system = new Ising2D(simInfo);
      break;

    default :
      std::cerr << "Physical system not specified. \n";
      std::cerr << "Aborting...\n";
      exit(10);
  }

}


#endif
