#include <iostream>
#include <fstream>
#include <string>             // std::string
#include <sstream>            // std::istringstream
#include "InputOutput.hpp"
//#include "RandomNumberGenerator.hpp"


// Input file should be specified as the first argument in command line
void readMainInputFile(const char* FileName, SimulationInfo& sim_info) {
  
  std::cout << "Reading main input file: " << FileName << std::endl;

  std::ifstream inputFile(FileName);
  std::string line, key;

  while (std::getline(inputFile, line)) {
    if (!line.empty()) {
      std::istringstream lineStream(line);
      lineStream >> key;

      if (key.compare(0, 1, "#") != 0) {

        // Read and set SimulationInfo
        if (key == "RestartFlag") {
          lineStream >> sim_info.restartFlag;
          std::cout << "Simulation Info: restartFlag = " << sim_info.restartFlag << std::endl;
          continue;
        }
        if (key == "PhysicalSystem") {
          lineStream >> sim_info.system;
          std::cout << "Simulation Info: system = " << sim_info.system << std::endl;
          continue;
        }
        if (key == "PhysicalSystemCommandLine") {
//          lineStream >> sim_info.physicalSystemCommandLine;
          lineStream.getline(sim_info.physicalSystemCommandLine, 256);
          std::cout << "Simulation Info: command line = " << sim_info.physicalSystemCommandLine << std::endl;
          continue;
        }
        if (key == "Algorithm") {
          lineStream >> sim_info.algorithm;
          std::cout << "Simulation Info: algorithm = " << sim_info.algorithm << std::endl;
          continue;
        }
        if (key == "RngSeed"){
          lineStream >> sim_info.rngSeed;
          std::cout << "Random number seed = " << sim_info.rngSeed << std::endl;
          continue;
        }
        if (key == "SpinModelLatticeSize") {
          lineStream >> sim_info.spinModelLatticeSize;
          std::cout << "Simulation Info: lattice size = " << sim_info.spinModelLatticeSize << std::endl;
          continue;
        }
        if (key == "QENumberOfAtoms") {
          lineStream >> sim_info.numAtoms;
          std::cout << "Simulation Info: number of atoms for Quantum Espresso = " << sim_info.numAtoms << std::endl;
          continue;
        }
        if (key == "MCInputFile") {
          lineStream >> sim_info.MCInputFile;
          std::cout << "Simulation Info: Input file for MC simulation = " << sim_info.MCInputFile << std::endl;
          continue;
        }
        if (key == "NumberOfRandomWalkers") {
          lineStream >> sim_info.numWalkers;
          std::cout << "Simulation Info: Number of random walkers = " << sim_info.numWalkers << std::endl;
          continue;
        }
        if (key == "NumberOfMPIranksPerWalker") {
          lineStream >> sim_info.numMPIranksPerWalker;
          std::cout << "Simulation Info: Number of MPI ranks random walker = Number of MPI ranks per physical system = " << sim_info.numMPIranksPerWalker << std::endl;
          continue;
        }

        // The following are inputs specific to different SimulationInfo parameters
        //switch (sim_info.system) {

        //  case 1 :
        //    std::cout << "Now read inputs for Quantum Espresso system." << std::endl;
        //    // specify QE input files
        //    break;

        //  case 2 :
        //    std::cout << "Now read inputs for LSMS system." << std::endl;
        //    // specify LSMS input files
        //    break;
 
        //  case 3 :
        //    std::cout << "Now read inputs for a Heisenberg 2D system." << std::endl;
        //    break;

        //  case 4 :
        //    std::cout << "Now read inputs for an Ising 2D system." << std::endl;
        //    break;

        //  default:
        //    {;}
        //}
        std::cout << "Unknown key: " << key << "  in input file " << FileName << std::endl;
        
      }
    }
  }

}

