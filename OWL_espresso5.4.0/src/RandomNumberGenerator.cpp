#include <iostream>
#include "RandomNumberGenerator.hpp"
#include "Communications.hpp"

std::mt19937 rng;
//int RngSeed {-1};

void initializeRandomNumberGenerator(int RngSeed)
{
 
  // get random seed
  if (RngSeed == -1) {
    RngSeed = int(time(NULL));
    std::cout << "Rank: " << myMPIRank
              << " No random number seed supplied. Take current time as a seed.\n";
  }
  else
    std::cout << "Random number seed supplied: " << RngSeed << std::endl;

  // YingWai: need to revise the communicator later when DFTCommunicator is created.
  MPI_Bcast(&RngSeed, 1, MPI_INT, 0, mpiCommunicator);
  std::cout << "Rank: " << myMPIRank << " Random number seed supplied: " << RngSeed << std::endl;

  // Initialize random number generator
  rng.seed(RngSeed);

}
