#include <iostream>
#include "RandomNumberGenerator.hpp"


// YingWai: should the RNG be placed as a member in MCAlgorithm class?
std::mt19937 rng;
//int RngSeed {-1};

void initializeRandomNumberGenerator(MPICommunicator phy_sys_comm, int RngSeed)
{
 
  // if random number seed is not supplied, get one from time 
  if (RngSeed == -1) {
    RngSeed = int(time(NULL));
    std::cout << "Rank: " << phy_sys_comm.thisMPIrank
              << " No random number seed supplied. Take current time as a seed.\n";
  }
  //else
    //std::cout << "Random number seed supplied: " << RngSeed << std::endl;

//YingWai's idea for optimization: the Bcast can be saved if every rank knows the group leader's ID
//                                 such that RngSeed += groupleader's ID

  // Make the seeds different for each walker
  if (phy_sys_comm.thisMPIrank == 0) {
    RngSeed += GlobalComm.thisMPIrank;
    std::cout << "Rank: " << GlobalComm.thisMPIrank << " Random number seed supplied: " << RngSeed << std::endl;
  }
  // Broadcast (synchronize) the random number seed within the same walker
  //MPI_Bcast(&RngSeed, 1, MPI_INT, 0, phy_sys_comm.communicator);
  phy_sys_comm.broadcastScalar(RngSeed, 0);

  // Initialize random number generators
  rng.seed(RngSeed);

}
