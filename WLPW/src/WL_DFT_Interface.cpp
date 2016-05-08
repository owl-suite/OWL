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
# include <mpi.h>
# include <iostream>
# include <fstream>
# include "MCAlgorithms.hpp"
using namespace std;


int main (int argc, char **argv) {
  
  // Set up MPI communicator:
  int exit_status;        // Environmental parameter
  int comm_help;          // Communicator
  MPI_Init(&argc, &argv);
  comm_help = MPI_Comm_c2f(MPI_COMM_WORLD);

  //YingWai's check   (Apr 7: can be removed when things work fine)
  YingWaisCheck(comm_help, exit_status);

  WangLandauSampling(comm_help, exit_status);

  return 0;
}
