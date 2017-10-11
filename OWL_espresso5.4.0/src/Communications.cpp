#include <iostream>
#include <vector>
#include "Communications.hpp"
#include "OWL_DFT_Interface.hpp"

// Global MPI communicator  (this is an extern variable)
MPICommunicator GlobalComm;


void MPICommunicator::initialize(MPI_Comm incomingComm)
{

  communicator = incomingComm;
  MPI_Comm_rank(communicator, &thisMPIrank);
  MPI_Comm_size(communicator, &totalMPIranks);
  //std::cout << "MPICommunicator: global rank = " << GlobalComm.thisMPIrank << ", local thisMPIrank = " << thisMPIrank << std::endl;

}


void MPICommunicator::finalize()
{
  // how to finalize just a Comm but not MPI_COMM_WORLD?
  // is it even necessary?
}


void MPICommunicator::barrier()
{

  MPI_Barrier(communicator);

}



void MPICommunicator::swapVector(void* data_ptr, int nElements, MPI_Datatype MPI_config_type, int partner)
{
 
  MPI_Status status;
  MPI_Sendrecv_replace(data_ptr, nElements, MPI_config_type, partner, 1, partner, 1, communicator, &status);

}



/// This function initializes the MPI communicators for global use, physical system, and the MC algorithm
void initializeMPICommunication(MPICommunicator& PhysicalSystemComm, 
                                MPICommunicator& MCAlgorithmComm) 
{

  MPI_Comm phys_sys_comm;
  MPI_Comm mc_alg_comm;

  /// 1. initialize global MPI communicator

  MPI_Init(NULL, NULL);
  GlobalComm.initialize(MPI_COMM_WORLD);

  /// 2. Define the MPI communicators for PhysicalSystem

  /// YingWai's note: Cannot use std::vector for walkerLeadersID because it does not match MPI_Group_incl's API
  //std::vector<int> walkerLeadersID;
  //walkerLeadersID.assign(simInfo.numWalkers, -1);
  int* walkerLeadersID;
  walkerLeadersID = new int[simInfo.numWalkers] {-1};

  /// exit if the total number of MPI ranks is not consistent with input info
  if (simInfo.numWalkers * simInfo.numMPIranksPerWalker != GlobalComm.totalMPIranks) { 

    std::cout << "ERROR!! Total number of MPI ranks is not consistent with input info.\n"
              << "        Please check the OWL input file or submission script.\n"
              << "        OWL aborting...\n";
    exit(7);

    /// ... or, assign the ranks here automatically
    /// TO DO: activate this when the basics are working fine  (Jun 24, 17)
    // simInfo.numMPIranksPerWalker = (GlobalComm.totalMPIranks - (GlobalComm.totalMPIranks % simInfo.numWalkers)) / simInfo.numWalkers; 

  }
  else {

    /// group processors into different walkers
    int walkerID = (GlobalComm.thisMPIrank - (GlobalComm.thisMPIrank % simInfo.numMPIranksPerWalker)) / simInfo.numMPIranksPerWalker;

    /// store the global MPI rank IDs for all "group leaders" (rank 0s) in each walker
    for (int i=0; i<simInfo.numWalkers; i++)
    {
      walkerLeadersID[i] = i * simInfo.numMPIranksPerWalker;
      if (GlobalComm.thisMPIrank == 0)
        std::cout << "walkerLeadersID[" << i << "] = " << walkerLeadersID[i] << "\n";
    }

    /// build the physical system communicators
    /// i.e., group the processors having the same walkerID
    MPI_Comm_split(MPI_COMM_WORLD, walkerID, 0, &phys_sys_comm);
    PhysicalSystemComm.initialize(phys_sys_comm);

    if (PhysicalSystemComm.thisMPIrank == 0)
      if (PhysicalSystemComm.totalMPIranks != simInfo.numMPIranksPerWalker) {
        std::cout << "ERROR!! Number of MPI ranks in a physical system communicator != number of MPI ranks per walker specified in OWL input file! \n "
                  << "        PhysicalSystemComm (walker) ID: " << walkerID << "\n"
                  << "        Number of MPI ranks assigned : " 
                  << PhysicalSystemComm.totalMPIranks << "\n"
                  << "        OWL aborting...\n";
        exit(7);
      }

  }

  
  /// 3. Define the MPI communicators for MCAlgorithm (this is the no-master (REWL) mode)
  MPI_Group MPI_GROUP_WORLD, WalkersGroup;
  MPI_Comm_group (MPI_COMM_WORLD, &MPI_GROUP_WORLD);
  MPI_Group_incl (MPI_GROUP_WORLD, simInfo.numWalkers, walkerLeadersID, &WalkersGroup);
  MPI_Comm_create (MPI_COMM_WORLD, WalkersGroup, &mc_alg_comm);
  if (PhysicalSystemComm.thisMPIrank == 0)
    MCAlgorithmComm.initialize(mc_alg_comm);

  delete walkerLeadersID;

}


void finalizeMPICommunication()
{

  MPI_Finalize();

}


// The following two should be moved to the Physical System MPI class, which in terms should be a member of the QuantumEspressoSystem class, which is then called in the QuantumEspressoSystem constructor
// i.e., initializeQEMPICommunication --> initializePhysicalSystemMPICommunication
//       finalizeQEMPICommunication   --> finalizePhysicalSystemMPICommunication
void initializeQEMPICommunication()
{
  //int comm_help = MPI_Comm_c2f(MPI_COMM_WORLD);   // MPI communicator handle for Fortran
  //owl_qe_startup(&comm_help);                     // Set up the PWscf calculation
  
  std::cout << "Intialized QE MPI communications..." << std::endl;
}


void finalizeQEMPICommunication()
{
  //int exit_status;                                // Environmental parameter for QE
  //owl_qe_stop(&exit_status);                      // Finish the PWscf calculation
  
  std::cout << "Finalized QE MPI communications..." << std::endl;
}

////////////////////////////////////////////////////////////////////////

