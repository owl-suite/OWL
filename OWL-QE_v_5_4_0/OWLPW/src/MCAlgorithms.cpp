#include <cstdio>
#include <cmath>
#include "MCAlgorithms.hpp"


void YingWaisCheck(int comm_help, int &exit_status)
{

  int natom;        // Total number of atoms in system
  double f_etot;    // Total energy of system (in Ry)
 
  initializeRandomNumberGenerator();
 
  wl_qe_startup_(&comm_help);        // Set up the PWscf calculation
  run_pwscf_(&exit_status);          // Execute the PWscf calculation

  //if (myMPIRank == 0) {
  get_natom_ener_(&natom, &f_etot);  // Extract the number of atoms and energy
 
  Matrix<double> pos_array;  // Set up the position array (in angstrom)
  pos_array.resize(3,natom); 
 
  Matrix<double> cell_array;  // Set up the cell vector array (in angstrom)
  cell_array.resize(3,3);
 
  get_pos_array_(&pos_array(0,0));  // Extract the position array from QE
  get_cell_array_(&cell_array(0,0));  // Extract the cell array from QE

  if (myMPIRank == 0) {
    // Write out system info
    writeSystemFile("energyFromMCAlgorithms.txt", f_etot, pos_array, cell_array);
  }
  proposeMCmoves(pos_array, cell_array);
  //MPI_Bcast(&pos_array(0,0), 3*natom, MPI_DOUBLE, 0, mpiCommunicator);
  //MPI_Bcast(&cell_array(0,0), 3*3, MPI_DOUBLE, 0, mpiCommunicator);
  //MPI_Barrier(mpiCommunicator);

  wl_stop_run_(&exit_status);  // Clean up the PWscf run

  // Based on the type of moves, one of these subroutines will be called 
  // to update either the atomic positions or lattice cell vector of  
  // simulation cell (needs to be merged with the Wang-Landau codes later...):

  pass_pos_array_(&pos_array(0,0));  // Update the atomic positions
  pass_cell_array_(&cell_array(0,0));  // Update the lattice cell vector

  // Call the Quantum Espresso subroutines to run the second PWscf calculation:
  wl_do_pwscf_(&exit_status);  // Run the subsequent PWscf calculation
  get_natom_ener_(&natom, &f_etot);

  get_pos_array_(&pos_array(0,0));
  get_cell_array_(&cell_array(0,0));

  //if (myMPIRank == 0)
    writeSystemFile("energy2.txt", f_etot, pos_array, cell_array);

  wl_qe_stop_(&exit_status);  // Finish the PWscf calculation
}



void WangLandauSampling(int comm_help, int &exit_status, int restartFlag)
{

  double startTime, currentTime, lastBackUpTime;
  startTime = lastBackUpTime = MPI_Wtime();
  if (myMPIRank == 0) 
    printf("WangLandauSampling is called.\n");

  // These should become Model class members 
  // System information
  int natom;                       // Total number of atoms in system
  double trialEnergy;              // Trial total energy of system (in Ry)
  double oldEnergy;                // Old total energy of the system (in Ry)
  Matrix<double> trialPos;         // Trial atom position    (in angstrom)
  Matrix<double> oldPos;           // Old atom position for restoration after a rejected MC move  (in angstrom)
  Matrix<double> trialLatticeVec;  // Trial cell vector   (in angstrom)
  Matrix<double> oldLatticeVec;    // Old cell vector for restoration after a rejected MC move (in angstrom)
  bool accepted = true;            // store the decision of WL acceptance for each move
  char fileName[51];
    
  // Initialize the simulation
  Histogram h(restartFlag);
  initializeRandomNumberGenerator();

  wl_qe_startup_(&comm_help);              // Set up the PWscf calculation
  run_pwscf_(&exit_status);                // Execute the PWscf calculation

  get_natom_ener_(&natom, &oldEnergy);     // Extract the number of atoms and energy
   
  trialPos.resize(3,natom);                // Resize position and cell vector arrays
  oldPos.resize(3,natom);
  trialLatticeVec.resize(3,3);
  oldLatticeVec.resize(3,3);
   
  get_pos_array_(&oldPos(0,0));            // Extract the position array from QE
  get_cell_array_(&oldLatticeVec(0,0));    // Extract the cell array from QE
  trialEnergy = oldEnergy;
  trialPos = oldPos;
  trialLatticeVec = oldLatticeVec;
 
  // Write out the energy
  if (myMPIRank == 0) 
    writeSystemFile("energyLatticePos.dat", oldEnergy, trialPos, trialLatticeVec);
    
  // Always accept the first energy
  h.updateHistogramDOS(oldEnergy);

  //YingWai: what is the purpose of this?
  wl_stop_run_(&exit_status);              // Clean up the PWscf run

//-------------- End initialization --------------//

// WL procedure starts here
  while (h.modFactor > h.modFactorFinal) {
    h.histogramFlat = false;

    while (!(h.histogramFlat)) {
      for (int MCSteps=0; MCSteps<h.histogramCheckInterval; MCSteps++)
      {

        proposeMCmoves(trialPos, trialLatticeVec);
        pass_pos_array_(&trialPos(0,0));          // Update the atomic positions
        pass_cell_array_(&trialLatticeVec(0,0));  // Update the lattice cell vector

        wl_do_pwscf_(&exit_status);               // Run the subsequent PWscf calculation 
        get_natom_ener_(&natom, &trialEnergy);

        // check if the energy falls within the energy range
        if ( !h.checkEnergyInRange(trialEnergy) )
          accepted = false;
        else {
          // determine WL acceptance
          if ( exp(h.getDOS(oldEnergy) - h.getDOS(trialEnergy)) > getRandomNumber() )
            accepted = true;
          else
            accepted = false;
        }

        if (accepted) {
           // Update histogram and DOS with trialEnergy
           h.updateHistogramDOS(trialEnergy);
           h.acceptedMoves++;        
 
           // Store trialPos, trialLatticeVec, trialEnergy
           oldPos = trialPos;
           oldLatticeVec = trialLatticeVec;
           oldEnergy = trialEnergy;
           if (myMPIRank == 0) {
             writeSystemFile("energyLatticePos.dat", oldEnergy, oldPos, oldLatticeVec);
             std::cerr << "trial move accepted: Energy = " << trialEnergy << "\n";
           }
        }
        else {
           if (myMPIRank == 0) {
             std::cerr << "trial move rejected: Energy = " << trialEnergy << "\n";
             std::cerr << "Update with energy = " << oldEnergy << "\n";
           }
           // Restore trialPos and trialLatticeVec
           trialPos = oldPos;
           trialLatticeVec = oldLatticeVec;

           // Update histogram and DOS with oldEnergy
           h.updateHistogramDOS(oldEnergy);
           h.rejectedMoves++;
        }
        h.totalMCsteps++;
      
        wl_stop_run_(&exit_status);              // Clean up the PWscf run

        // Write restart files at interval
        currentTime = MPI_Wtime();
        if (myMPIRank == 0) {
          if (currentTime - lastBackUpTime > 300) {
            h.writeHistogramDOSFile("hist_dos_checkpoint.dat");
            writeQErestartFile("OWL_QE_restart_input", trialPos, trialLatticeVec);
            lastBackUpTime = currentTime;
          }
        }
      }

      // Check histogram flatness
      h.histogramFlat = h.checkHistogramFlatness();
      //currentTime = MPI_Wtime();
      //if (currentTime - lastBackUpTime > 1650) {
      //  h.writeHistogramDOSFile("hist_dos_checkpoint.dat");
      //  writeQErestartFile("OWL_QE_restart_input", trialPos, trialLatticeVec);
      //  lastBackUpTime = currentTime;
      //}
    }

    // Go to next iteration
    h.modFactor /= h.modFactorReducer;
    h.resetHistogram();
    h.iterations++;

    if (myMPIRank == 0) {
      printf("Number of iterations performed = %d\n", h.iterations);
    
      // Also write restart file here 
      sprintf(fileName, "hist_dos_iteration%02d.dat", h.iterations);
      h.writeHistogramDOSFile(fileName);
      writeQErestartFile("OWL_QE_restart_input", trialPos, trialLatticeVec);
    }
  }

  wl_qe_stop_(&exit_status);  // Finish the PWscf calculation

}
