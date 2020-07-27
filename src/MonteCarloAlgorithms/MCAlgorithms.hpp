#ifndef MC_ALGORITHMS_HPP
#define MC_ALGORITHMS_HPP

#include <limits>
#include <mpi.h>
#include "PhysicalSystems/PhysicalSystemBase.hpp"
#include "Main/Globals.hpp"


enum OutputMode {checkPoint, endOfIteration, endOfSimulation};

// Base class for all Monte Carlo algorithms
class MonteCarloAlgorithm {

public :

  // Constructor
  MonteCarloAlgorithm();

  // Destructor
  virtual ~MonteCarloAlgorithm() {}

  virtual void run() = 0;

  // MPI Communicator for communications among MC walkers
  //MPICommunicator MCAlgorithmCommunicator;

protected :

  double currentTime;
  double lastBackUpTime;
  double checkPointInterval {600.0};                           // unit: seconds.

  unsigned long int configurationWriteInterval { std::numeric_limits<unsigned long int>::max() };

  // it stores the decision of acceptance for each move
  bool acceptMove {false};

  // MC statistics:
  unsigned long int totalMCsteps  {0};
  unsigned long int acceptedMoves {0};
  unsigned long int rejectedMoves {0};

private :


};

#endif
