#include <iostream>
#include "Communications.hpp"
#include "OWL_DFT_Interface.hpp"

//Global variables for MPI info (these are the extern variables)
int myMPIRank;
int numMPIRanks;
MPI_Comm mpiCommunicator;
//MPI_Status mpiStatus;
//int exit_status;
//int comm_help;


void initializeMPICommunication(SimulationInfo simInfo) 
{
  MPI_Init(NULL, NULL);

  mpiCommunicator = MPI_COMM_WORLD;
  MPI_Comm_rank(mpiCommunicator, &myMPIRank);
  MPI_Comm_size(mpiCommunicator, &numMPIRanks);
  std::cout << "myMPIrank = " << myMPIRank << std::endl;

  // Here, the MPI communicators for (i) MCAlgorithm and (ii) PhysicalSystem should be defined

/*
  switch (simInfo.system) {
    case 1 :
      //initializeQEMPICommunication();
      //comm_help = MPI_Comm_c2f(MPI_COMM_WORLD);
      //owl_qe_startup(&comm_help);        // Set up the PWscf calculation
      break;

    case 2 :
      // initializeLSMSMPICommunication();
      break;

    default :  // this includes cases 3 and 4 (all simple models)
      // initializeSimpleMPICommunication();
      // basically nothing; 1 MPI per system
      {};
  }
*/

}


void finalizeMPICommunication(SimulationInfo simInfo)
{

/*
  switch (simInfo.system) {
    case 1 :
      finalizeQEMPICommunication();
      //owl_qe_stop(&exit_status);  // Finish the PWscf calculation
      break;

    case 2 :
      //finalizeLSMSMPICommunication();
      break;

    default :  // this includes cases 3 and 4 (all simple models)
      MPI_Finalize();   // this might need to be moved outside the switch

  }
*/

  MPI_Finalize();

}

// The following two should be moved to the Physical System MPI class, which in terms should be a member of the QuantumEspressoSystem class, which is then called in the QuantumEspressoSystem constructor
// i.e., initializeQEMPICommunication --> initializePhysicalSystemMPICommunication
//       finalizeQEMPICommunication   --> finalizePhysicalSystemMPICommunication
void initializeQEMPICommunication()
{
  // !!! MPI_COMM_WORLD should be changed to the MPI Comm. group after the MC Algorithm Communicator is set up. !!!
  //int comm_help = MPI_Comm_c2f(MPI_COMM_WORLD);   // MPI communicator handle for Fortran
  //owl_qe_startup(&comm_help);                     // Set up the PWscf calculation
  
  std::cout << "Intialized QE MPI communications..." << std::endl;
  std::cout << "myMPIrank = " << myMPIRank << std::endl;
}


void finalizeQEMPICommunication()
{
  //int exit_status;                                // Environmental parameter for QE
  //owl_qe_stop(&exit_status);                      // Finish the PWscf calculation
  
  std::cout << "Finalized QE MPI communications..." << std::endl;
  std::cout << "myMPIrank = " << myMPIRank << std::endl;
}
