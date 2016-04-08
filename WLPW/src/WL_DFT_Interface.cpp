//
// This file is distributed under the terms of the
// GNU General Public License. See the file `License'
// in the root directory of the present distribution,
// or http://www.gnu.org/copyleft/gpl.txt.
//
//-----------------------------------------------------------------
//   WL+DFT Interface Program
//-----------------------------------------------------------------
//
// Main C++ program calling the Quantum Espresso subroutines for
// performing PWscf calculations and passing the energy, position,
// and cell vector arrays to the Wang-Landau algorithm.
//
# include <mpi.h>
# include <iostream>
# include <fstream>
# include "MCAlgorithms.hpp"
# include "WL_DFT_Interface.hpp"

using namespace std;


int main (int argc, char **argv) {
  int natom;        // A total number of atoms in system
  double f_etot;    // Total energy of system (in eV/unit cell)
  int exit_status;  // Environmental parameter
  int comm_help;    // Communicator

  // Set up the communicator
  MPI_Init(&argc, &argv);
  comm_help = MPI_Comm_c2f(MPI_COMM_WORLD);

  wl_qe_startup_(&comm_help);        // Set up the PWscf calculation
  run_pwscf_(&exit_status);          // Execute the PWscf calculation
  get_natom_ener_(&natom, &f_etot);  // Extract the number of atoms and energy

  //YingWai's check   (Apr 7: can be removed when things work fine)
  YingWaisCheck();

  // Write out the energy
  writeEnergyFile("energy.txt", f_etot);

  Matrix<double> pos_array;  // Set up the position array (in angstrom)
  pos_array.resize(3,natom);

  Matrix<double> cell_array;  // Set up the cell vector array (in angstrom)
  cell_array.resize(3,3);

  get_array_(&pos_array(0,0), &cell_array(0,0));

  proposeMCmoves(pos_array, cell_array);

  wl_stop_run_(&exit_status);  // Clean up the PWscf run
  wl_pass_array_(&pos_array(0,0));  // Pass the position array to QE

  wl_do_pwscf_(&exit_status);  // Run the subsequent PWscf calculation

  get_natom_ener_(&natom, &f_etot);

  writeEnergyFile("energy2.txt", f_etot);

  get_array_(&pos_array(0,0), &cell_array(0,0));

  wl_qe_stop_(&exit_status);  // Finish the PWscf calculation
}
