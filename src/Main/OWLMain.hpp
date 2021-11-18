#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include "Globals.hpp"
#include "MonteCarloAlgorithms/MCAlgorithms.hpp"
#include "MonteCarloAlgorithms/Metropolis.hpp"
#include "MonteCarloAlgorithms/WangLandauSampling.hpp"
#include "MonteCarloAlgorithms/ReplicaExchangeWangLandau.hpp"
#include "MonteCarloAlgorithms/MulticanonicalSampling.hpp"
#include "MonteCarloAlgorithms/HistogramFreeMUCA.hpp"
#include "PhysicalSystems/Heisenberg2D.hpp"
#include "PhysicalSystems/Heisenberg3D.hpp"
#include "PhysicalSystems/Ising2D.hpp"
#include "PhysicalSystems/CrystalStructure3D.hpp"
#include "PhysicalSystems/Alloy3D.hpp"
#include "PhysicalSystems/HeisenbergHexagonal2D.hpp"

#ifdef DRIVER_MODE_QE
#include "PhysicalSystems/QuantumEspresso/QuantumEspressoSystem.hpp"
#endif

void readCommandLineArguments(int argc, char* argv[]) {

  // Read in input file information
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " [OWL's Input File Name] \n";
    exit(7);
  }

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
      std::cerr << "In standalone mode, use of Quantum Espresso is not supported.\n"
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

    case 5 :
      physical_system = new Heisenberg3D();
      break;

    case 6 :
      physical_system = new CrystalStructure3D(simInfo.MCInputFile);
      break;

    case 7 :
      physical_system = new Alloy3D(simInfo.MCInputFile);
      break;

  case 8 :
    physical_system = new HeisenbergHexagonal2D();
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
      MC = new Metropolis( physical_system, simInfo.MCInputFile);
      break;

    case 2 :
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
      MC = new WangLandauSampling( physical_system );
  }

}


void finalizeSimulation(PhysicalSystem*      &physical_system,
                        MonteCarloAlgorithm* &MC)
{
  
  delete physical_system;
  delete MC;

}


#endif
