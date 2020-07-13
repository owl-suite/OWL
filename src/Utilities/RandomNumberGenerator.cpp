#include <iostream>
#include "RandomNumberGenerator.hpp"


// YingWai: should the RNG be placed as a member in MCAlgorithm class?
std::mt19937 rng_engine;
//int RngSeed {-1};
std::uniform_real_distribution<double> distribution1(-0.5, 0.5); 
std::uniform_real_distribution<double> distribution2(0.0, 1.0); 
std::uniform_int_distribution<int> distribution_int;
std::uniform_int_distribution<unsigned int> distribution_uint; 


void initializeRandomNumberGenerator(MPICommunicator phy_sys_comm, int RngSeed)
{
 
  // if random number seed is not supplied, get one from time 
  if (RngSeed == -1) {
    RngSeed = int(time(NULL));
    std::cout << "   Random number seed       :  " << RngSeed << "   (from current system time)\n";
  }
  else
    std::cout << "   Random number seed       :  " << RngSeed << "   (from input file)\n\n";

//YingWai's idea for optimization: the Bcast can be saved if every rank knows the group leader's ID
//                                 such that RngSeed += groupleader's ID

  // Make the seeds different for each walker
  if (phy_sys_comm.thisMPIrank == 0) {
    RngSeed += GlobalComm.thisMPIrank;
    //std::cout << "Rank: " << GlobalComm.thisMPIrank << " Random number seed supplied: " << RngSeed << std::endl;
  }
  // Broadcast (synchronize) the random number seed within the same walker
  //MPI_Bcast(&RngSeed, 1, MPI_INT, 0, phy_sys_comm.communicator);
  phy_sys_comm.broadcastScalar(RngSeed, 0);

  // Initialize random number generators
  rng_engine.seed(RngSeed);

}
