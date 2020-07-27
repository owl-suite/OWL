#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <string>             // std::string
#include <sstream>            // std::istringstream
#include "InputOutput.hpp"


// Input file should be specified as the first argument in command line
void readMainInputFile(const char* FileName)
{
  
  //std::cout << "Reading main input file: " << FileName << std::endl;
  int length = strlen(FileName);
  simInfo.MCInputFile = new char[length+1]();
  strncpy(simInfo.MCInputFile, FileName, length+1);
  //std::cout << "simInfo.MCInputFile set to be: " << simInfo.MCInputFile << std::endl;

  std::ifstream inputFile(FileName);
  std::string line, key;

  if (inputFile.is_open()) {

    while (std::getline(inputFile, line)) {
      if (!line.empty()) {
        std::istringstream lineStream(line);
        lineStream >> key;
  
        if (key.compare(0, 1, "#") != 0) {
  
          // Read and set SimulationInfo
          if (key == "RestartFlag") {
            lineStream >> simInfo.restartFlag;
            //std::cout << "Simulation Info: restartFlag = " << simInfo.restartFlag << std::endl;
            continue;
          }
          else if (key == "PhysicalSystem") {
            lineStream >> simInfo.system;
            //std::cout << "Simulation Info: system = " << simInfo.system << std::endl;
            continue;
          }
          else if (key == "PhysicalSystemCommandLine") {
  //          lineStream >> simInfo.physicalSystemCommandLine;
            lineStream.getline(simInfo.physicalSystemCommandLine, 256);
            //std::cout << "Simulation Info: command line = " << simInfo.physicalSystemCommandLine << std::endl;
            continue;
          }
          else if (key == "Algorithm") {
            lineStream >> simInfo.algorithm;
            //std::cout << "Simulation Info: algorithm = " << simInfo.algorithm << std::endl;
            continue;
          }
          else if (key == "RngSeed"){
            lineStream >> simInfo.rngSeed;
            //std::cout << "Random number seed = " << simInfo.rngSeed << std::endl;
            continue;
          }
          else if (key == "SpinModelLatticeSize") {
            lineStream >> simInfo.spinModelLatticeSize;
            //std::cout << "Simulation Info: lattice size = " << simInfo.spinModelLatticeSize << std::endl;
            continue;
          }
          else if (key == "QENumberOfAtoms") {
            lineStream >> simInfo.numAtoms;
            //std::cout << "Simulation Info: number of atoms for Quantum Espresso = " << simInfo.numAtoms << std::endl;
            continue;
          }
          else if (key == "QEMCMoveSet") {
            lineStream >> simInfo.QEMCMoveSet;
            //std::cout << "Simulation Info: choice of Monte Carlo move set for Quantum Espresso system= " << simInfo.QEMCMoveSet << std::endl;
            continue;
          }
          else if (key == "numberOfWindows") {
            lineStream >> simInfo.numberOfWindows;
            //std::cout << "Simulation Info: number of energy window in Replica-Exchange Wang-Landau = " << simInfo.numberOfWindows << std::endl;
            continue;
          }
          else if (key == "numberOfWalkersPerWindow") {
            lineStream >> simInfo.numberOfWalkersPerWindow;
            //std::cout << "Simulation Info: number of random walkers per window in Replica-Exchange Wang-Landau = " << simInfo.numberOfWalkersPerWindow << std::endl;
            continue;
          }
          else if (key == "NumberOfMPIranksPerWalker") {
            lineStream >> simInfo.numMPIranksPerWalker;
            //std::cout << "Simulation Info: Number of MPI ranks random walker = Number of MPI ranks per physical system = " << simInfo.numMPIranksPerWalker << std::endl;
            continue;
          }
        
        }

      }
    }
 
    inputFile.close();
  }

  if ((simInfo.numberOfWindows > 0) && (simInfo.numberOfWalkersPerWindow > 0))
    simInfo.numWalkers = simInfo.numberOfWindows * simInfo.numberOfWalkersPerWindow;
 
}


void writeSimulationInfo()
{

  printf("\n");
  printf("                  ####################################################################\n");
  printf("                  ##                                                                ##\n");
  printf("                  ##   Oak-Ridge/Open-source Wang-Landau (OWL) simulation package   ##\n");
  printf("                  ##   Copyright (C) 2015-2020                                      ##\n");
  printf("                  ##                                                                ##\n");
  printf("                  ####################################################################\n");
  printf("\n\n");

  printf("Initializing simulation...\n");


  switch (simInfo.restartFlag) {
    case 0 : 
      // fresh start
      std::cout << "   New/restarted simulation :  New\n";
      break;
    case 1 : 
      // restart from file
      std::cout << "   New/restarted simulation :  Restarted\n";
      break;
    default :
      std::cout << "   New/restarted simulation :  Restart flag != 0 or 1. Use default: new simulation.\n";
      simInfo.restartFlag = 0;
    }


  switch (simInfo.system) {
    case 1 :
      std::cout << "   Physical system          :  Quantum Espresso (driver mode)\n";
      std::cout << "    - Command line options:    " << simInfo.physicalSystemCommandLine << "\n";
      break;

    case 2 :
      std::cout << "   Physical system          :  Locally Self-Consistent Multiple Scattering (driver mode)\n";
      std::cout << "    - Command line options:    " << simInfo.physicalSystemCommandLine << "\n";
      break;

    case 3 :
      std::cout << "   Physical system          :  2D Heisenberg model\n";
      break;

    case 4 :
      std::cout << "   Physical system          :  2D Ising model\n";
      break;

    case 5 :
      std::cout << "   Physical system          :  3D Heisenberg model\n";
      break;

    case 6 :
      std::cout << "   Physical system          :  3D Customized crystal structure\n";
      break;

    default :
      std::cerr << "   Physical system          :  ERROR! Physical system not specified. \n";
      std::cerr << "\nOWL Aborting...\n";
      exit(10);
  }


  switch (simInfo.algorithm) {
    case 1 :
      std::cout << "   MC algorithm             :  Metropolis Sampling\n";
      break;

    case 2 :
      std::cout << "   MC algorithm             :  Wang-Landau Sampling\n";
      break;

    case 3 :
      std::cout << "   MC algorithm             :  Multicanonical Sampling\n";
      break;

    case 4 :
      std::cout << "   MC algorithm             :  Parallel Tempering\n";
      break;

    case 5 :
      std::cout << "   MC algorithm             :  Replica-Exchange Wang-Landau sampling\n";
      break;
      
    case 6 :
      std::cout << "   MC algorithm             :  Discrete Histogram-Free Multicanonical Sampling\n";
      break;
      
    default :
      std::cout << "   MC algorithm             :  MC algorithm not specified. Use default: Wang-Landau sampling.\n";
  }

}