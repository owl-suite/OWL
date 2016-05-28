#include "InputOutput.hpp"

void writeEnergyFile(char fileName[], double energy)
{ 
  FILE *energy_file;
  energy_file = fopen(fileName, "a");
  fprintf(energy_file, "%f", energy);
  fclose(energy_file); 
} 

void writeQErestartFile(char fileName[], Matrix<double> atom_positions)
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
  fprintf(QE_file, "    ibrav= 6\n");
  fprintf(QE_file, "    celldm(1) = 7.28330256917748\n");
  fprintf(QE_file, "    celldm(3) = 1.04661740543287\n");
  fprintf(QE_file, "    nat= 5\n");
  fprintf(QE_file, "    ntyp= 3\n");
  fprintf(QE_file, "    ecutwfc = 50\n");
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
  fprintf(QE_file, "ATOMIC_POSITIONS {crystal}\n");
  fprintf(QE_file, "Pb    %14.9f %14.9f %14.9f\n", atom_positions(0,0), 
                                                   atom_positions(1,0), 
                                                   atom_positions(2,0) );
  fprintf(QE_file, "Ti    %14.9f %14.9f %14.9f\n", atom_positions(0,1),
                                                   atom_positions(1,1), 
                                                   atom_positions(2,1) );
  fprintf(QE_file, "O     %14.9f %14.9f %14.9f\n", atom_positions(0,2),
                                                   atom_positions(1,2), 
                                                   atom_positions(2,2) );
  fprintf(QE_file, "O     %14.9f %14.9f %14.9f\n", atom_positions(0,3),
                                                   atom_positions(1,3), 
                                                   atom_positions(2,3) );
  fprintf(QE_file, "O     %14.9f %14.9f %14.9f\n", atom_positions(0,4),
                                                   atom_positions(1,4), 
                                                   atom_positions(2,4) );
  fprintf(QE_file, "\n");
  fprintf(QE_file, "K_POINTS automatic\n");
  fprintf(QE_file, "4 4 4 0 0 0\n");

  fclose(QE_file);
}
