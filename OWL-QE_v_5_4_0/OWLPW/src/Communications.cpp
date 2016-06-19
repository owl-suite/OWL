#include <iostream>
#include "Communications.hpp"

int myMPIRank;
int numMPIRanks;
MPI_Comm mpiCommunicator;
MPI_Status mpiStatus;

void initializeMPICommunication() 
{
  MPI_Init(NULL, NULL);
  mpiCommunicator = MPI_COMM_WORLD;
  MPI_Comm_rank(mpiCommunicator, &myMPIRank);
  MPI_Comm_size(mpiCommunicator, &numMPIRanks);
 
  std::cout << "myMPIrank = " << myMPIRank << std::endl;
}

void finalizeMPICommunication()
{
  MPI_Finalize();
}
