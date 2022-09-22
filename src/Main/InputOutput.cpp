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
  
  std::cout << "Reading main input file: " << FileName << "\n\n";
  
  size_t length = strlen(FileName);
  simInfo.MCInputFile = new char[length+1]();
  strncpy(simInfo.MCInputFile, FileName, length+1);
  //std::cout << "simInfo.MCInputFile set to be: " << simInfo.MCInputFile << "\n";

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
            //std::cout << "Simulation Info: restartFlag = " << simInfo.restartFlag << "\n";
            continue;
          }
          else if (key == "PhysicalSystem") {
            lineStream >> simInfo.system;
            //std::cout << "Simulation Info: system = " << simInfo.system << "\n";
            continue;
          }
          else if (key == "PhysicalSystemCommandLine") {
  //          lineStream >> simInfo.physicalSystemCommandLine;
            lineStream.getline(simInfo.physicalSystemCommandLine, 256);
            //std::cout << "Simulation Info: command line = " << simInfo.physicalSystemCommandLine << "\n";
            continue;
          }
          else if (key == "Algorithm") {
            lineStream >> simInfo.algorithm;
            //std::cout << "Simulation Info: algorithm = " << simInfo.algorithm << "\n";
            continue;
          }
          else if (key == "RngSeed"){
            lineStream >> simInfo.rngSeed;
            //std::cout << "Random number seed = " << simInfo.rngSeed << "\n";
            continue;
          }
          else if (key == "SpinModelLatticeSize") {
            lineStream >> simInfo.spinModelLatticeSize;
            //std::cout << "Simulation Info: lattice size = " << simInfo.spinModelLatticeSize << "\n";
            continue;
          }
          else if (key == "SpinConfigInitMethod") {
            lineStream >> simInfo.spinConfigInitMethod;
            //std::cout << "Simulation Info: method to initialize spin configuration = " << simInfo.spinConfigInitMethod << "\n";
            continue;
          }
          else if (key == "SpinModelDimension") {
            lineStream >> simInfo.spinModelDimension;
            //std::cout << "Simulation Info: lattice size = " << simInfo.spinModelLatticeSize << "\n";
            continue;
          }
          else if (key == "QENumberOfAtoms") {
            lineStream >> simInfo.numAtoms;
            //std::cout << "Simulation Info: number of atoms for Quantum Espresso = " << simInfo.numAtoms << "\n";
            continue;
          }
          else if (key == "QEMCMoveSet") {
            lineStream >> simInfo.QEMCMoveSet;
            //std::cout << "Simulation Info: choice of Monte Carlo move set for Quantum Espresso system= " << simInfo.QEMCMoveSet << "\n";
            continue;
          }
          else if (key == "numberOfWindows") {
            lineStream >> simInfo.numberOfWindows;
            //std::cout << "Simulation Info: number of energy window in Replica-Exchange Wang-Landau = " << simInfo.numberOfWindows << "\n";
            continue;
          }
          else if (key == "numberOfWalkersPerWindow") {
            lineStream >> simInfo.numberOfWalkersPerWindow;
            //std::cout << "Simulation Info: number of random walkers per window in Replica-Exchange Wang-Landau = " << simInfo.numberOfWalkersPerWindow << "\n";
            continue;
          }
          else if (key == "NumberOfMPIranksPerWalker") {
            lineStream >> simInfo.numMPIranksPerWalker;
            //std::cout << "Simulation Info: Number of MPI ranks random walker = Number of MPI ranks per physical system = " << simInfo.numMPIranksPerWalker << "\n";
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


void writeBanner()
{

  printf("\n");
  printf("                  ####################################################################\n");
  printf("                  ##                                                                ##\n");
  printf("                  ##   Open-source/Oak-Ridge Wang-Landau (OWL) simulation package   ##\n");
  printf("                  ##   Copyright (C) 2015-2020                                      ##\n");
  printf("                  ##                                                                ##\n");
  printf("                  ####################################################################\n");
  printf("\n\n");

}


void writeSimulationInfo()
{

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

    case 7 :
      std::cout << "   Physical system          :  3D Alloy system\n";
      break;

    case 8 :
      std::cout << "   Physical system          :  2D Hexagonal Heisenberg model\n";
      break; 

    case 9 :
      std::cout << "   Physical system          :  ND Ising model\n";
      break; 

    case 10 :
      std::cout << "   Physical system          :  2D Ising model with next nearest neighbor interactions\n";
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
