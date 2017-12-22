#include <iostream>
#include <fstream>
#include <cstring>
#include <string>             // std::string
#include <sstream>            // std::istringstream
#include "InputOutput.hpp"


// Input file should be specified as the first argument in command line
void readMainInputFile(const char* FileName, SimulationInfo& sim_info) {
  
  std::cout << "Reading main input file: " << FileName << std::endl;
  int length = strlen(FileName);
  sim_info.MCInputFile = new char[length +1]();
  strncpy(sim_info.MCInputFile, FileName, length); 
  //std::cout << "sim_info.MCInputFile set to be: " << sim_info.MCInputFile << std::endl;

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
            lineStream >> sim_info.restartFlag;
            //std::cout << "Simulation Info: restartFlag = " << sim_info.restartFlag << std::endl;
            continue;
          }
          else if (key == "PhysicalSystem") {
            lineStream >> sim_info.system;
            //std::cout << "Simulation Info: system = " << sim_info.system << std::endl;
            continue;
          }
          else if (key == "PhysicalSystemCommandLine") {
  //          lineStream >> sim_info.physicalSystemCommandLine;
            lineStream.getline(sim_info.physicalSystemCommandLine, 256);
            //std::cout << "Simulation Info: command line = " << sim_info.physicalSystemCommandLine << std::endl;
            continue;
          }
          else if (key == "Algorithm") {
            lineStream >> sim_info.algorithm;
            //std::cout << "Simulation Info: algorithm = " << sim_info.algorithm << std::endl;
            continue;
          }
          else if (key == "RngSeed"){
            lineStream >> sim_info.rngSeed;
            //std::cout << "Random number seed = " << sim_info.rngSeed << std::endl;
            continue;
          }
          else if (key == "SpinModelLatticeSize") {
            lineStream >> sim_info.spinModelLatticeSize;
            //std::cout << "Simulation Info: lattice size = " << sim_info.spinModelLatticeSize << std::endl;
            continue;
          }
          else if (key == "QENumberOfAtoms") {
            lineStream >> sim_info.numAtoms;
            //std::cout << "Simulation Info: number of atoms for Quantum Espresso = " << sim_info.numAtoms << std::endl;
            continue;
          }
          else if (key == "NumberOfRandomWalkers") {
            lineStream >> sim_info.numWalkers;
            //std::cout << "Simulation Info: Number of random walkers = " << sim_info.numWalkers << std::endl;
            continue;
          }
          else if (key == "NumberOfMPIranksPerWalker") {
            lineStream >> simInfo.numMPIranksPerWalker;
            //std::cout << "Simulation Info: Number of MPI ranks random walker = Number of MPI ranks per physical system = " << sim_info.numMPIranksPerWalker << std::endl;
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
}

