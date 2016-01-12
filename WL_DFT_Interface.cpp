//
// This file is distributed under the terms of the
// GNU General Public License. See the file `License'
// in the root directory of the present distribution,
// or http://www.gnu.org/copyleft/gpl.txt.
//
//--------------------------------------------------------------
// WL+DFT Interface Program
//--------------------------------------------------------------
//
// Main C++ program to call the appropriate subroutines from
// Quantum Espresso codes for performing PWscf calculations and
// pass the energy and position arrays to Wang-Landau algorithm.
//
# include <iostream>
# include <fstream>
# include "Matrix.hpp"  // Matrix class header

using namespace std;

// Define the external Quantum Espresso F90 subroutines

extern "C" {
 void run_qe_startup_( );
 void run_qe_run_(int *exit_status);
 void run_qe_stop_(int *exit_status);
 void get_natom_ener_(int *natom, double *f_etot);
 void get_pos_array_(double *pos_array);
 void get_cell_array_(double *cell_array);
}

// Call the Quantum Espresso subroutines to run PWscf calculations.

int main () {
  int natom;        // A total number of atoms in system
  double f_etot;    // Total energy of system (in eV/unit cell)
  int exit_status;  // Environmental parameter

  run_qe_startup_( );
  run_qe_run_(&exit_status);
  get_natom_ener_(&natom, &f_etot);  // Call and copy the total number of atoms and energy

  FILE *n_atom_file;
  n_atom_file = fopen("n_atom.txt", "w");  // Write out the total number of atoms
  fprintf(n_atom_file, "%i", natom);
  fclose(n_atom_file);

  FILE *energy_file;
  energy_file = fopen("energy.txt", "w");  // Write out the energy
  fprintf(energy_file, "%f", f_etot);
  fclose(energy_file);

  Matrix<double> pos_array;  // Set up the position array (in angstrom)
  pos_array.resize(3,natom);
  get_pos_array_(&pos_array(0,0));

  FILE *pos_file;
  pos_file = fopen("pos_array.txt", "w");  // Write out the position array
  for (int i = 0; i < natom; i++)
    {fprintf(pos_file, "%f %f %f \n", pos_array(0,i), pos_array(1,i), pos_array(2,i));
    }
  fclose(pos_file);

  Matrix<double> cell_array;  // Set up the cell vector array (in angstrom)
  cell_array.resize(3,3);
  get_cell_array_(&cell_array(0,0));

  FILE *cell_file;
  cell_file = fopen("cell_array.txt", "w");  // Write out the cell vector array
  for (int j = 0; j < 3; j++)
    {fprintf(cell_file, "%f %f %f \n", cell_array(0,j), cell_array(1,j), cell_array(2,j));
    }
  fclose(cell_file);

  run_qe_stop_(&exit_status);  // Finish the PWscf calculations and close all the files
}
