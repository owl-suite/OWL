#ifndef COMMUNICATIONS_HPP
#define COMMUNICATIONS_HPP

#include "mpi.h"
#include "Globals.hpp"

// Global variables for MPI info
extern int myMPIRank;                 // MPI rank for this processor
extern int numMPIRanks;               // Total number of processors
extern MPI_Comm mpiCommunicator;      // Global Communicator
//extern MPI_Status mpiStatus;        // MPI status for the communicator
//extern int exit_status;               // Environmental parameter for QE
//extern int comm_help;                 // MPI communicator handle for Fortran


void initializeMPICommunication(SimulationInfo);
void finalizeMPICommunication(SimulationInfo);

void initializeQEMPICommunication();
void finalizeQEMPICommunication();

/////////////////////////////////////////////////////////////////////////////////////////
// Ying Wai's plan: MPI communications should be rearranged with the following classes //
/////////////////////////////////////////////////////////////////////////////////////////

/*
class MPICommunicatorBase {

public :

  MPI_Comm comm;
  int myRank;
  int totalRanks;

  void initializeCommunicator();
  void finalizeCommunicator();

};



class MCAlgorithmCommunicator : public MPICommunicatorBase 
{

};


class PhysicalSystemCommunicator : public MPICommunicatorBase
{

};
*/


#endif
