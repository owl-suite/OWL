#include <iostream>
#include "Communications.hpp"
#include "WL_DFT_Interface.hpp"

//Global variables for MPI info (these are the extern variables)
int myMPIRank;
int numMPIRanks;
MPI_Comm mpiCommunicator;
//MPI_Status mpiStatus;
//int exit_status;
//int comm_help;


void initializeMPICommunication() 
{
  //MPI_Init(NULL, NULL);
  //comm_help = MPI_Comm_c2f(MPI_COMM_WORLD);
  //wl_qe_startup_(&comm_help);        // Set up the PWscf calculation

  mpiCommunicator = MPI_COMM_WORLD;
  MPI_Comm_rank(mpiCommunicator, &myMPIRank);
  MPI_Comm_size(mpiCommunicator, &numMPIRanks);
  std::cout << "myMPIrank = " << myMPIRank << std::endl;
}

void finalizeMPICommunication()
{
  //wl_qe_stop_(&exit_status);  // Finish the PWscf calculation
  //MPI_Finalize();
}


