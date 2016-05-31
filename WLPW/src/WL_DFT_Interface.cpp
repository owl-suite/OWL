// Markus: Prefer Modified BSD license. CHECK!!  (May 13, 16)
//         Dual license? 
//
// This file is distributed under the terms of the
// GNU General Public License. See the file `License'
// in the root directory of the present distribution,
// or http://www.gnu.org/copyleft/gpl.txt.
//
//--------------------------------------------------------------------
//   WL+DFT Interface Program
//--------------------------------------------------------------------
//
// Main C++ program calling the Quantum Espresso subroutines for
// performing PWscf calculations and passing the energy, position,
// and cell vector arrays to the Wang-Landau algorithm.
//
# include <iostream>
# include <fstream>
# include "MCAlgorithms.hpp"

int main (int argc, char *argv[]) {

  // Set up MPI communicator:
  int exit_status;        // Environmental parameter
  int comm_help;          // Communicator
  MPI_Init(&argc, &argv);
  comm_help = MPI_Comm_c2f(MPI_COMM_WORLD);

  // Read in restarting information
  int restartFlag {-1};
  if (argc != 2)
    std::cout << "Usage: " << argv[0] << " [RestartFlag=0/1] \n";
  else {
    if (sscanf(argv[1], "%d", &restartFlag) != 1) {
      std::cout << "Error: cannot read restart flag! \n";
      exit(1);
    }
    else
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
  }

  //YingWai's check   (Apr 7: can be removed when things work fine)
  //YingWaisCheck(comm_help, exit_status);

  WangLandauSampling(comm_help, exit_status, restartFlag);




  return 0;
}
