#include <cstring>
#include <cstdio>
#include "QuantumEspressoSystem.hpp"
#include "QEMCMoves.hpp"
#include "OWL_DFT_Interface.hpp"
#include "Main/Communications.hpp"

//Constructor
QuantumEspressoSystem::QuantumEspressoSystem(MPICommunicator PhysicalSystemComm)
{
  natom       = simInfo.numAtoms;  // TO DO: this should be cross-checked with QE input file!
  //oldEnergy   = 0.0;
  //trialEnergy = 0.0;
  
  initializeObservables(1); // observable[0] = energy

  readCommandLineOptions();

  //initializeQEMPICommunication();
  int comm_help = MPI_Comm_c2f(PhysicalSystemComm.communicator);  // MPI communicator handle for Fortran
  owl_qe_startup(&comm_help, &nimage, &npool, &ntg, &nband, &ndiag, QEInputFile);  // Set up the PWscf calculation
  
  if (GlobalComm.thisMPIrank == 0)
    std::cout << "Initialized QE MPI communications..." << std::endl;

  run_pwscf_(&MPI_exit_status);                                // Execute the PWscf calculation
  get_natom_ener(&natom, &observables[0]);                     // Extract the number of atoms and energy
  //observables[0] = trialEnergy;
  //std::cout << "Here: " << natom << ", " << trialEnergy << std::endl;

  trialConfig.atomic_positions.resize(3,natom);                // Resize positions, cell vectors and atomic species
  trialConfig.lattice_vectors.resize(3,3);
  trialConfig.atomic_species.resize(natom);
  oldConfig.atomic_positions.resize(3,natom);
  oldConfig.lattice_vectors.resize(3,3);
  oldConfig.atomic_species.resize(natom);

  buildMPIConfigurationType();
  pointerToConfiguration = static_cast<void*>(&trialConfig);   // Is it ok to point to a struct like this? (Dec 26, 17)

  // check if the following should be called in here:
  // yes, because they initialize trialConfig.atomic_positions and trialConfig.lattice_vectors
  get_pos_array(&trialConfig.atomic_positions(0,0));           // Extract the position array from QE
  get_cell_array(&trialConfig.lattice_vectors(0,0));           // Extract the cell array from QE
  get_atomic_species(&trialConfig.atomic_species[0]);          // Exrract the atomic species from QE
  owl_stop_run(&MPI_exit_status);                              // Clean up the PWscf run

}

// Destructor
QuantumEspressoSystem::~QuantumEspressoSystem()
{

  // finalizeQEMPICommunication();
  int exit_status;                                // Environmental parameter for QE
  owl_qe_stop(&exit_status);                      // Finish the PWscf calculation

  // Free MPI datatype
  pointerToConfiguration = NULL;
  MPI_Type_free(&MPI_ConfigurationType);

  std::cout << "Finalized QE MPI communications...\n";

  deleteObservables();

  std::cout << "QuantumEspressoSystem finished\n";

}


void QuantumEspressoSystem::readCommandLineOptions()
{
  
  // TODO: this should be PhysicalSystemComm after Issue #7 is taken care of.  (Feb 13, 19)
  if (GlobalComm.thisMPIrank == 0) {
    std::cout << "Reading the following command line for Quantum Espresso: \"" 
              << simInfo.physicalSystemCommandLine << "\"" << std::endl;
  }

  char* pch;
  pch = strtok (simInfo.physicalSystemCommandLine, " ");
  while (pch != NULL)
  {
    //printf ("%s\n", pch);

        if (strncmp("-i",pch,2) == 0) {
            pch = strtok (NULL, " ");
            //printf ("%s\n", pch);
            strncpy(QEInputFile, pch, 80);
            QEInputFile[80] = '\0';
            pch = strtok (NULL, " ");
            continue;
        }   
        if ((strncmp("-ni",pch,3) == 0) 
             || (strncmp("-npot",pch,5) == 0)) {
            pch = strtok (NULL, " ");
            //printf ("%s\n", pch);
            nimage = std::atoi(pch);
            pch = strtok (NULL, " ");
            continue;
        }   
        //if (strncmp("-npot",pch,5) == 0) {
        //    pch = strtok (NULL, " ");
        //    printf ("%s\n", pch);
        //    npots = std::atoi(pch);
        //    pch = strtok (NULL, " ");
        //    continue;
        //}   
        if ((strncmp("-nk",pch,3) == 0)
             || (strncmp("-npoo",pch,5) == 0)) {
            pch = strtok (NULL, " ");
            //printf ("%s\n", pch);
            npool = std::atoi(pch);
            pch = strtok (NULL, " ");
            continue;
        }   
        if (strncmp("-nt",pch,3) == 0) {
            pch = strtok (NULL, " ");
            //printf ("%s\n", pch);
            ntg = std::atoi(pch);
            pch = strtok (NULL, " ");
            continue;
        }
        if (strncmp("-nb",pch,3) == 0) {
            pch = strtok (NULL, " ");
            //printf ("%s\n", pch);
            nband = std::atoi(pch);
            pch = strtok (NULL, " ");
            continue;
        }
        if ((strncmp("-nd",pch,3) == 0)
            || (strncmp("-no",pch,3) == 0)
            || (strcmp("-nproc_diag",pch) == 0)
            || (strcmp("-nproc_ortho",pch) == 0)) {
            pch = strtok (NULL, " ");
            //printf ("%s\n", pch);
            ndiag = std::atoi(pch);
            pch = strtok (NULL, " ");
            continue;
        }
        //if (strncmp("-nr",pch,3) == 0) {
        //    pch = strtok (NULL, " ");
        //    printf ("%s\n", pch);
        //    nres = std::atoi(pch);
        //    pch = strtok (NULL, " ");
        //    continue;
        //}
        std::cerr << "Error when reading QE command line options!!\n "
                  << std::endl;

  }  

};


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

  owl_do_pwscf(&MPI_exit_status);               // Run the subsequent PWscf calculation
  get_natom_ener(&natom, &observables[0]);      // Obtain the # of atoms and energy from QE
  //observables[0] = trialEnergy;
  owl_stop_run(&MPI_exit_status);               // Clean up the PWscf run

}


void QuantumEspressoSystem::doMCMove()
{

  proposeMCmoves(trialConfig.atomic_positions, trialConfig.lattice_vectors, trialConfig.atomic_species);
  pass_pos_array(&trialConfig.atomic_positions(0,0));       // Update the atomic positions
  pass_cell_array(&trialConfig.lattice_vectors(0,0));       // Update the lattice cell vector
  pass_atomic_species(&trialConfig.atomic_species[0]);      // Update the atomic species
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
  //oldEnergy = trialEnergy;

  // update "old" configurations
  oldConfig.lattice_vectors  = trialConfig.lattice_vectors;
  oldConfig.atomic_positions = trialConfig.atomic_positions;
  oldConfig.atomic_species   = trialConfig.atomic_species;
 
}



void QuantumEspressoSystem::rejectMCMove()
{

  // Restore trialConfig.atomic_positions and trialConfig.lattice_vectors
  for (int i=0; i<numObservables; i++)
    observables[i] = oldObservables[i];
  //trialEnergy = oldEnergy;
  trialConfig.lattice_vectors  = oldConfig.lattice_vectors;
  trialConfig.atomic_positions = oldConfig.atomic_positions;
  trialConfig.atomic_species   = oldConfig.atomic_species;

}
 

void QuantumEspressoSystem::writeSystemFile(const char* fileName)
{
    FILE *system_file;
    system_file = fopen(fileName, "a");
    fprintf(system_file, " %14.9f ", observables[0]);
 
     for (unsigned int i=0; i<oldConfig.lattice_vectors.n_col(); i++)
       for(unsigned int j=0; j<oldConfig.lattice_vectors.n_row(); j++)
         fprintf(system_file, " %14.9f ", oldConfig.lattice_vectors(j,i));
 
     for (unsigned int i=0; i<oldConfig.atomic_positions.n_col(); i++)
       for(unsigned int j=0; j<oldConfig.atomic_positions.n_row(); j++)
         fprintf(system_file, " %14.9f ", oldConfig.atomic_positions(j,i));

     for(unsigned int i=0; i<oldConfig.atomic_species.size(); i++)
       fprintf(system_file, " %3d ", oldConfig.atomic_species[i]);

     fprintf(system_file, "\n");
     fclose(system_file);
}


void QuantumEspressoSystem::writeQErestartFile(const char* fileName)
{
  
  // Now hard-coded, need to generalize in the future.
  // Should look if QE has this functionality already.

/*
  FILE *QE_file;
  QE_file = fopen(fileName, "w");

  fprintf(QE_file, "&control\n");
  fprintf(QE_file, "   calculation = 'scf'\n");
  fprintf(QE_file, "   restart_mode = 'from_scratch'\n");
  fprintf(QE_file, "   forc_conv_thr = 1.0d-4\n");
  fprintf(QE_file, "   forc_conv_thr = 3.0d-4\n");
  fprintf(QE_file, "   tstress = .true.\n");
  fprintf(QE_file, "   tprnfor = .true.\n");
  fprintf(QE_file, "   pseudo_dir = './'\n");
  fprintf(QE_file, "!   lberry = .true.\n");
  fprintf(QE_file, "!   gdir = 3\n");
  fprintf(QE_file, "!   nppstr = 8 \n");
  fprintf(QE_file, "/\n");
  fprintf(QE_file, "&system\n");
  fprintf(QE_file, "    ibrav = 0\n");
  //fprintf(QE_file, "!   celldm(1) = 1.0\n");
  fprintf(QE_file, "    nat= 5\n");
  fprintf(QE_file, "    ntyp= 3\n");
  fprintf(QE_file, "    ecutwfc = 50\n");
  fprintf(QE_file, "    ecutrho = 400\n");
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
  fprintf(QE_file, "    cell_factor = 30.0d0\n");
  fprintf(QE_file, "/\n");
  fprintf(QE_file, "\n");
  fprintf(QE_file, "ATOMIC_SPECIES\n");
  fprintf(QE_file, "Pb  207.2    Pb.pz-d-van.UPF\n");
  fprintf(QE_file, "Ti  47.867   022-Ti-ca-sp-vgrp_serge.uspp\n");
  fprintf(QE_file, "O   16.00    O_ps.uspp.UPF\n");
  fprintf(QE_file, "\n");
  fprintf(QE_file, "ATOMIC_POSITIONS {angstrom}\n");
  fprintf(QE_file, "Pb    %14.9f %14.9f %14.9f\n", oldConfig.atomic_positions(0,0), 
                                                   oldConfig.atomic_positions(1,0), 
                                                   oldConfig.atomic_positions(2,0) );
  fprintf(QE_file, "Ti    %14.9f %14.9f %14.9f\n", oldConfig.atomic_positions(0,1),
                                                   oldConfig.atomic_positions(1,1), 
                                                   oldConfig.atomic_positions(2,1) );
  fprintf(QE_file, "O     %14.9f %14.9f %14.9f\n", oldConfig.atomic_positions(0,2),
                                                   oldConfig.atomic_positions(1,2), 
                                                   oldConfig.atomic_positions(2,2) );
  fprintf(QE_file, "O     %14.9f %14.9f %14.9f\n", oldConfig.atomic_positions(0,3),
                                                   oldConfig.atomic_positions(1,3), 
                                                   oldConfig.atomic_positions(2,3) );
  fprintf(QE_file, "O     %14.9f %14.9f %14.9f\n", oldConfig.atomic_positions(0,4),
                                                   oldConfig.atomic_positions(1,4), 
                                                   oldConfig.atomic_positions(2,4) );
  fprintf(QE_file, "\n");
  fprintf(QE_file, "K_POINTS automatic\n");
  fprintf(QE_file, "8 8 8 0 0 0\n");
  fprintf(QE_file, "\n");
  fprintf(QE_file, "CELL_PARAMETERS {angstrom}\n");
  for (unsigned int i=0; i<oldConfig.lattice_vectors.n_col(); i++) {
    for(unsigned int j=0; j<oldConfig.lattice_vectors.n_row(); j++)
      fprintf(QE_file, " %14.9f ", oldConfig.lattice_vectors(j,i));
    fprintf(QE_file, "\n");
  }

  fclose(QE_file);

*/

}


void QuantumEspressoSystem::buildMPIConfigurationType()
{

  int          count {3};
  int          numElements[count];
  MPI_Aint     address_displacements[count];
  MPI_Datatype types[count];

  MPI_Aint start_address;
  MPI_Aint address1;
  MPI_Aint address2;

  // Number of elements in each array
  numElements[0] = trialConfig.atomic_positions.size();
  numElements[1] = trialConfig.lattice_vectors.size();
  numElements[2] = trialConfig.atomic_species.size();
  //std::cout << "Debugging check: numElements[0]/[1] = " << numElements[0] << " , " << numElements[1] << std::endl;

  // The derived datatype consists of two arrays of MPI_DOUBLE
  types[0] = types[1] = MPI_DOUBLE;
  types[2] = MPI_INT;

  // Calculate the relative addresses for each array
  MPI_Get_address(&trialConfig.atomic_positions[0], &start_address);
  MPI_Get_address(&trialConfig.lattice_vectors[0], &address1);
  MPI_Get_address(&trialConfig.atomic_species[0], &address2);
  address_displacements[0] = 0;
  address_displacements[1] = address1 - start_address;
  address_displacements[2] = address2 - start_address;

  MPI_Type_create_struct(count, numElements, address_displacements, types, &MPI_ConfigurationType);
  MPI_Type_commit(&MPI_ConfigurationType);

}

