#include "mpi.h"

// Global variables for MPI info
extern int myMPIRank;                 // MPI rank for this processor
extern int numMPIRanks;               // Total number of processors
extern MPI_Comm mpiCommunicator;      // Global Communicator
//extern MPI_Status mpiStatus;        // MPI status for the communicator
//extern int exit_status;               // Environmental parameter for QE
//extern int comm_help;                 // MPI communicator handle for Fortran


void initializeMPICommunication();
void finalizeMPICommunication();

