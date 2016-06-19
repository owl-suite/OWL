#include "InputOutput.hpp"

void writeSystemFile(char fileName[], double energy, Matrix<double> atom_positions, Matrix<double> cell_vectors)
{ 
  FILE *system_file;
  system_file = fopen(fileName, "a");
  fprintf(system_file, " %14.9f ", energy);

  for (int i=0; i<cell_vectors.n_col(); i++)
    for(int j=0; j<cell_vectors.n_row(); j++)
      fprintf(system_file, " %14.9f ", cell_vectors(j,i));
  
  for (int i=0; i<cell_vectors.n_col(); i++)
    for(int j=0; j<cell_vectors.n_row(); j++)
      fprintf(system_file, " %14.9f ", atom_positions(j,i));

  fprintf(system_file, "\n");
  fclose(system_file); 
} 

void writeQErestartFile(char fileName[], Matrix<double> atom_positions, 
                        Matrix<double> cell_vectors)
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
  fprintf(QE_file, "\n");
  fprintf(QE_file, "CELL_PARAMETERS {angstrom}\n");
  for (int i=0; i<cell_vectors.n_col(); i++) {
    for(int j=0; j<cell_vectors.n_row(); j++)
      fprintf(QE_file, " %14.9f ", cell_vectors(j,i));
    fprintf(QE_file, "\n");
  }

  fclose(QE_file);
}
