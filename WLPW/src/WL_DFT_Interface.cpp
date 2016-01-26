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
# include <iostream>
# include <fstream>
# include "Matrix.hpp"   // Matrix class header

using namespace std;

// Define the external Quantum Espresso F90 subroutines

extern "C" {
  void wl_qe_startup_( );
  void run_pwscf_(int *exit_status);
  void wl_qe_stop_(int *exit_status);
  void get_natom_ener_(int *natom, double *f_etot);
  void get_array_(double *pos_array, double *cell_array);
  void wl_pass_array_(double *pos_array);
  void wl_do_pwscf_(int *exit_status);
  void wl_stop_run_(int *exit_status);
}

// Call the Quantum Espresso subroutines to run PWscf calculations.

int main () {
  int natom;        // A total number of atoms in system
  double f_etot;    // Total energy of system (in eV/unit cell)
  int exit_status;  // Environmental parameter

  wl_qe_startup_( );
  run_pwscf_(&exit_status);
  get_natom_ener_(&natom, &f_etot);  // Extract the number of atoms and energy

  FILE *energy_file;
  energy_file = fopen("energy.txt", "w");  // Write out the energy
  fprintf(energy_file, "%f", f_etot);
  fclose(energy_file);

  Matrix<double> pos_array;  // Set up the position array (in angstrom)
  pos_array.resize(3,natom);

  Matrix<double> cell_array;  // Set up the cell vector array (in angstrom)
  cell_array.resize(3,3);

  get_array_(&pos_array(0,0), &cell_array(0,0));

  wl_stop_run_(&exit_status);  // Clean up the PWscf run
  wl_pass_array_(&pos_array(0,0));  // Pass the position array to QE

  wl_do_pwscf_(&exit_status);  // Run the subsequent PWscf run

  get_natom_ener_(&natom, &f_etot);

  FILE *energy_file_2;
  energy_file_2 = fopen("energy_2.txt", "w");
  fprintf(energy_file_2, "%f", f_etot);
  fclose(energy_file_2);

  get_array_(&pos_array(0,0), &cell_array(0,0));

  wl_qe_stop_(&exit_status);  // Finish the PWscf calculations
}
