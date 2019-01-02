#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include "Globals.hpp"
#include "MCAlgorithms.hpp"
#include "WangLandauSampling.hpp"
#include "ReplicaExchangeWangLandau.hpp"
#include "MulticanonicalSampling.hpp"
#include "HistogramFreeMUCA.hpp"
#include "Heisenberg2D.hpp"
#include "Ising2D.hpp"

#ifdef DRIVER_MODE_QE
#include "QuantumEspresso/QuantumEspressoSystem.hpp"
#endif

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

void setSimulation(PhysicalSystem*      &physical_system,
                   MonteCarloAlgorithm* &MC,
                   MPICommunicator      physicalSystemComm,
                   MPICommunicator      mcAlgorithmComm)
{
  
  // Determine Physical System 
  //  1: QuantumExpresso
  //  2: LSMS  
  //  3: Heisenberg 2D
  //  4: Ising 2D
  //  5: ...
  switch (simInfo.system) {
    case 1 :
#ifdef DRIVER_MODE_QE
      physical_system = new QuantumEspressoSystem( physicalSystemComm );
#else
      std::cerr << "In standalone mode, using Quantum Espresso is not supported.\n"
                << "Please recompile OWL with \'make owl-qe\'.\n";
#endif
      break;

    case 2 :
      std::cout << "LSMS not yet available.\n";
      break;

    case 3 :
      physical_system = new Heisenberg2D();
      break;

    case 4 :
      physical_system = new Ising2D();
      break;

    default :
      std::cerr << "Physical system not specified. \n";
      std::cerr << "Aborting...\n";
      exit(10);
  }

  // Determine MC algorithm
  //  1. Metropolis sampling
  //  2. Wang-Landau sampling
  //  3. Multicanonical sampling (MUCA)
  //  4. Parallel tempering
  //  5. Replica-Exchange Wang-Landau sampling (REWL)
  //  6. Histogram-free Multicanonical sampling (discrete energy version)
  switch (simInfo.algorithm) {
    case 1 :
      MC = new Metropolis( physical_system );
      break;

    case 2 :
      //MC = new WangLandauSampling( simInfo.restartFlag, simInfo.MCInputFile );
      MC = new WangLandauSampling( physical_system );
      break;

    case 3 :
      MC = new MulticanonicalSampling( physical_system );
      break;

    case 4 :
      std::cout << "Parallel tempering not yet implemented\n";
      exit(10);
      break;

    case 5 :
      MC = new ReplicaExchangeWangLandau( physical_system, physicalSystemComm, mcAlgorithmComm );
      break;
      
    case 6 :
      MC = new DiscreteHistogramFreeMUCA( physical_system );
      break;
      
    default :
      std::cout << "Monte Carlo algorithm not specified. Use default: Wang-Landau sampling.\n";
  }

}


void finalizeSimulation(PhysicalSystem*      &physical_system,
                        MonteCarloAlgorithm* &MC)
{
  
  delete physical_system;
  delete MC;

}


#endif
