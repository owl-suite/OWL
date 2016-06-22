#include "mpi.h"

extern int myMPIRank;                 // MPI rank for this processor
extern int numMPIRanks;               // Total number of processors
extern MPI_Comm mpiCommunicator;      // Global Communicator
extern MPI_Status mpiStatus;          // MPI status for the communicator

void initializeMPICommunication();
void finalizeMPICommunication();

