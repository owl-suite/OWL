#include "QuantumEspressoSystem.hpp"
#include "MCMoves.hpp"
#include "Communications.hpp"
#include "OWL_DFT_Interface.hpp"

//Constructor
QuantumEspressoSystem::QuantumEspressoSystem(SimulationInfo& sim_info)
{
  natom       = sim_info.size;
  oldEnergy   = 0.0;
  trialEnergy = 0.0;
  
  initializeObservables(1); // observable[0] = energy

  std::cout << "Initializing Quantum Espresso with the following command line: \"" << sim_info.commandline << "\"" << std::endl;
  command_line_options_mp_get_command_line_(sim_info.commandline,256);
  //__command_line_options_MOD_get_command_line(sim_info.commandline,256);

//  int comm_help = MPI_Comm_c2f(MPI_COMM_WORLD);   // MPI communicator handle for Fortran
//  owl_qe_startup_(&comm_help);                    // Set up the PWscf calculation
//
//  std::cout << "Intialized QE MPI communications..." << std::endl;
//  std::cout << "myMPIrank = " << myMPIRank << std::endl;

  run_pwscf_(&MPI_exit_status);                 // Execute the PWscf calculation
  get_natom_ener_(&natom, &trialEnergy);        // Extract the number of atoms and energy
  std::cout << "Here: " << natom << ", " << trialEnergy << std::endl;
  observables[0] = trialEnergy;

  trialPos.resize(3,natom);                     // Resize position and cell vector arrays
  oldPos.resize(3,natom);
  trialLatticeVec.resize(3,3);
  oldLatticeVec.resize(3,3);

  //check if the following should be called in here:
  // yes, because they intialize trialPos and trialLatticeVec
  get_pos_array_(&trialPos(0,0));               // Extract the position array from QE
  get_cell_array_(&trialLatticeVec(0,0));       // Extract the cell array from QE
  owl_stop_run_(&MPI_exit_status);              // Clean up the PWscf run

}

//Destructor
QuantumEspressoSystem::~QuantumEspressoSystem()
{
//  int exit_status;                              // Environmental parameter for QE
//  owl_qe_stop_(&exit_status);                   // Finish the PWscf calculation

  std::cout << "Finalized QE MPI communications..." << std::endl;
  std::cout << "myMPIrank = " << myMPIRank << std::endl;

  deleteObservables();
}


void QuantumEspressoSystem::writeConfiguration(int option, const char* fileName)
{

  switch (option) {

  case 1 : { 
    writeQErestartFile(fileName);
    break;
  }
  default : {
    writeSystemFile(fileName);
  }

  }

}


void QuantumEspressoSystem::getObservables()
{
  // trialEnergy should be changed to observables[0]

  owl_do_pwscf_(&MPI_exit_status);               // Run the subsequent PWscf calculation
  get_natom_ener_(&natom, &trialEnergy);
  observables[0] = trialEnergy;
  owl_stop_run_(&MPI_exit_status);               // Clean up the PWscf run

}


void QuantumEspressoSystem::doMCMove()
{

  proposeMCmoves(trialPos, trialLatticeVec);
  pass_pos_array_(&trialPos(0,0));              // Update the atomic positions
  pass_cell_array_(&trialLatticeVec(0,0));      // Update the lattice cell vector
  
}

/*
void QuantumEspressoSystem::undoMCMove()
{

}
*/

void QuantumEspressoSystem::acceptMCMove()
{
  // update "old" observables
  for (int i=0; i<numObservables; i++)
    oldObservables[i] = observables[i];
    oldEnergy      = trialEnergy;

  // update "old" configurations
  oldLatticeVec  = trialLatticeVec;
  oldPos         = trialPos;
  
  // Old implementation of OWL-QE
  //trialEnergy     = oldEnergy;
  //trialPos        = oldPos;
  //trialLatticeVec = oldLatticeVec;
}



void QuantumEspressoSystem::rejectMCMove()
{

  // Restore trialPos and trialLatticeVec
  for (int i=0; i<numObservables; i++)
    observables[i] = oldObservables[i];
  trialEnergy      = oldEnergy;
  trialLatticeVec  = oldLatticeVec;
  trialPos         = oldPos;

}
 

void QuantumEspressoSystem::writeSystemFile(const char* fileName)
{
    FILE *system_file;
    system_file = fopen(fileName, "a");
    fprintf(system_file, " %14.9f ", observables[0]);
 
     for (unsigned int i=0; i<oldLatticeVec.n_col(); i++)
       for(unsigned int j=0; j<oldLatticeVec.n_row(); j++)
         fprintf(system_file, " %14.9f ", oldLatticeVec(j,i));
 
     for (unsigned int i=0; i<oldPos.n_col(); i++)
       for(unsigned int j=0; j<oldPos.n_row(); j++)
         fprintf(system_file, " %14.9f ", oldPos(j,i));
 
     fprintf(system_file, "\n");
     fclose(system_file);
}


void QuantumEspressoSystem::writeQErestartFile(const char* fileName)
{
  
  // Now hard-coded, need to generalize in the future.
  // Should look if QE has this functionality already.

  FILE *QE_file;
  QE_file = fopen(fileName, "w");

  fprintf(QE_file, "&control\n");
  fprintf(QE_file, "   calculation = 'scf'\n");
  fprintf(QE_file, "   restart_mode = 'from_scratch'\n");
  fprintf(QE_file, "   forc_conv_thr = 3.0d-4\n");
  fprintf(QE_file, "   tstress = .true.\n");
  fprintf(QE_file, "   tprnfor = .true.\n");
  fprintf(QE_file, "   pseudo_dir = './'\n");
  fprintf(QE_file, "!   lberry = .true.\n");
  fprintf(QE_file, "!   gdir = 3\n");
  fprintf(QE_file, "!   nppstr = 8 \n");
  fprintf(QE_file, "/\n");
  fprintf(QE_file, "&system\n");
  fprintf(QE_file, "    ibrav= 0\n");
  fprintf(QE_file, "!   celldm(1) = 1.0\n");
  fprintf(QE_file, "    nat= 5\n");
  fprintf(QE_file, "    ntyp= 3\n");
  fprintf(QE_file, "    ecutwfc = 50\n");
  fprintf(QE_file, "    nosym = .true.\n");
  fprintf(QE_file, "!    nspin = 2     ! 1 = non-polarized 2 = spin-polarized\n");
  fprintf(QE_file, "!    occupations = 'smearing'\n");
  fprintf(QE_file, "!    smearing = 'methfessel-paxton'\n");
  fprintf(QE_file, "!    degauss = 0.02\n");
  fprintf(QE_file, "!    starting_magnetization(1) = 2\n");
  fprintf(QE_file, "!    starting_magnetization(2) = -2\n");
  fprintf(QE_file, "/\n");
  fprintf(QE_file, "&electrons\n");
  fprintf(QE_file, "    electron_maxstep = 1000\n");
  fprintf(QE_file, "    mixing_beta = 0.7\n");
  fprintf(QE_file, "    conv_thr = 1.0d-8\n");
  fprintf(QE_file, "/\n");
  fprintf(QE_file, "&ions\n");
  fprintf(QE_file, "/\n");
  fprintf(QE_file, "&cell\n");
  fprintf(QE_file, "    cell_factor = 3.0d0\n");
  fprintf(QE_file, "/\n");
  fprintf(QE_file, "\n");
  fprintf(QE_file, "ATOMIC_SPECIES\n");
  fprintf(QE_file, "Pb  207.2    Pb.pz-d-van.UPF\n");
  fprintf(QE_file, "Ti  47.867   022-Ti-ca-sp-vgrp_serge.uspp\n");
  fprintf(QE_file, "O   16.00    O_ps.uspp.UPF\n");
  fprintf(QE_file, "\n");
  fprintf(QE_file, "ATOMIC_POSITIONS {angstrom}\n");
  fprintf(QE_file, "Pb    %14.9f %14.9f %14.9f\n", oldPos(0,0), 
                                                   oldPos(1,0), 
                                                   oldPos(2,0) );
  fprintf(QE_file, "Ti    %14.9f %14.9f %14.9f\n", oldPos(0,1),
                                                   oldPos(1,1), 
                                                   oldPos(2,1) );
  fprintf(QE_file, "O     %14.9f %14.9f %14.9f\n", oldPos(0,2),
                                                   oldPos(1,2), 
                                                   oldPos(2,2) );
  fprintf(QE_file, "O     %14.9f %14.9f %14.9f\n", oldPos(0,3),
                                                   oldPos(1,3), 
                                                   oldPos(2,3) );
  fprintf(QE_file, "O     %14.9f %14.9f %14.9f\n", oldPos(0,4),
                                                   oldPos(1,4), 
                                                   oldPos(2,4) );
  fprintf(QE_file, "\n");
  fprintf(QE_file, "K_POINTS automatic\n");
  fprintf(QE_file, "4 4 4 0 0 0\n");
  fprintf(QE_file, "\n");
  fprintf(QE_file, "CELL_PARAMETERS {angstrom}\n");
  for (unsigned int i=0; i<oldLatticeVec.n_col(); i++) {
    for(unsigned int j=0; j<oldLatticeVec.n_row(); j++)
      fprintf(QE_file, " %14.9f ", oldLatticeVec(j,i));
    fprintf(QE_file, "\n");
  }

  fclose(QE_file);

}
