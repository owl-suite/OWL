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
# include <mpi.h>
# include <iostream>
# include <fstream>
# include "MCAlgorithms.hpp"

int main (int argc, char **argv) {
  
  // Set up MPI communicator:
  int exit_status;        // Environmental parameter
  int comm_help;          // Communicator
  MPI_Init(&argc, &argv);
  comm_help = MPI_Comm_c2f(MPI_COMM_WORLD);
  //wl_qe_startup_(&comm_help);              // Set up the PWscf calculation

  //YingWai's check   (Apr 7: can be removed when things work fine)
  //YingWaisCheck(comm_help, exit_status);

  WangLandauSampling(comm_help, exit_status);




  //wl_qe_stop_(&exit_status);               // Finish the PWscf calculation
  return 0;
}
