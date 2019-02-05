#include <iostream>
#include <fstream>
#include <cstring>
#include <string>             // std::string
#include <sstream>            // std::istringstream
#include "InputOutput.hpp"


// Input file should be specified as the first argument in command line
void readMainInputFile(const char* FileName) {
  
  std::cout << "Reading main input file: " << FileName << std::endl;
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
          //else if (key == "NumberOfRandomWalkers") {
          //  lineStream >> simInfo.numWalkers;
            //std::cout << "Simulation Info: Number of random walkers = " << simInfo.numWalkers << std::endl;
          //  continue;
          //}
          else if (key == "numberOfWindows") {
            lineStream >> simInfo.numberOfWindows;
            continue;
          }
          else if (key == "numberOfWalkersPerWindow") {
            lineStream >> simInfo.numberOfWalkersPerWindow;
            continue;
          }
          else if (key == "NumberOfMPIranksPerWalker") {
            lineStream >> simInfo.numMPIranksPerWalker;
            //std::cout << "Simulation Info: Number of MPI ranks random walker = Number of MPI ranks per physical system = " << simInfo.numMPIranksPerWalker << std::endl;
            continue;
          }
          else {
            //std::cout << "Unknown key: " << key << "  in input file " << FileName << std::endl;
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

