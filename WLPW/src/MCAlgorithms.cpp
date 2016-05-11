#include "MCAlgorithms.hpp"


void YingWaisCheck(int comm_help, int exit_status)
{

  int natom;        // Total number of atoms in system
  double f_etot;    // Total energy of system (in Ry)
 
  initializeRandomNumberGenerator();
 
  wl_qe_startup_(&comm_help);        // Set up the PWscf calculation
  run_pwscf_(&exit_status);          // Execute the PWscf calculation
  get_natom_ener_(&natom, &f_etot);  // Extract the number of atoms and energy
 
  // Write out the energy
  writeEnergyFile("energyFromMCAlgorithms.txt", f_etot);

  Matrix<double> pos_array;  // Set up the position array (in angstrom)
  pos_array.resize(3,natom); 
 
  Matrix<double> cell_array;  // Set up the cell vector array (in angstrom)
  cell_array.resize(3,3);
 
  get_pos_array_(&pos_array(0,0));  // Extract the position array from QE
  get_cell_array_(&cell_array(0,0));  // Extract the cell array from QE

  proposeMCmoves(pos_array, cell_array);

  wl_stop_run_(&exit_status);  // Clean up the PWscf run

  // Based on the type of moves, one of these subroutines will be called 
  // to update either the atomic positions or lattice cell vector of  
  // simulation cell (needs to be merged with the Wang-Landau codes later...):

  pass_pos_array_(&pos_array(0,0));  // Update the atomic positions
  pass_cell_array_(&cell_array(0,0));  // Update the lattice cell vector

  // Call the Quantum Espresso subroutines to run the second PWscf calculation:
  wl_do_pwscf_(&exit_status);  // Run the subsequent PWscf calculation
  get_natom_ener_(&natom, &f_etot);

  writeEnergyFile("energy2.txt", f_etot);

  get_pos_array_(&pos_array(0,0));
  get_cell_array_(&cell_array(0,0));

  wl_qe_stop_(&exit_status);  // Finish the PWscf calculation
}



void WangLandauSampling(int comm_help, int exit_status)
{

  // These should become Model class members 
  // System information
  int natom;        // Total number of atoms in system
  double f_etot;    // Total energy of system (in Ry)


// Initialize the simulation
  Histogram h;
  h.getNumberOfBins();    // check, can be removed later

  initializeRandomNumberGenerator();
 
  wl_qe_startup_(&comm_help);        // Set up the PWscf calculation
  run_pwscf_(&exit_status);          // Execute the PWscf calculation
  get_natom_ener_(&natom, &f_etot);  // Extract the number of atoms and energy
 
  // Write out the energy
  writeEnergyFile("energyFromMCAlgorithms.txt", f_etot);

  Matrix<double> pos_array;  // Set up the position array (in angstrom)
  pos_array.resize(3,natom); 
 
  Matrix<double> cell_array;  // Set up the cell vector array (in angstrom)
  cell_array.resize(3,3);
 
  get_pos_array_(&pos_array(0,0));    // Extract the position array from QE
  get_cell_array_(&cell_array(0,0));  // Extract the cell array from QE

  
// WL procedure starts here


}
