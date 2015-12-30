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
// passing the energy and position arrays to Wang-Landau algorithm.
//
# include <iostream>

extern "C" {
 void getenergy_startup_( );
 void getenergy_run_(int *exit_status);
}

// Call the QE f90 subroutines to run PWscf calculations.

void main () {
  int exit_status;
  getenergy_startup_( );
  getenergy_run_(&exit_status);
}
