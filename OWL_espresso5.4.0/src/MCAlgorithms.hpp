#ifndef MC_ALGORITHMS_HPP
#define MC_ALGORITHMS_HPP

#include <mpi.h>
#include "PhysicalSystemBase.hpp"
#include "Globals.hpp"


// Base class for all Monte Carlo algorithms
class MonteCarloAlgorithm {

public :

  // Constructor
  MonteCarloAlgorithm();

  // Destructor
  virtual ~MonteCarloAlgorithm() {}

  virtual void run(PhysicalSystem*) = 0;

  // Should they be set directly by I/O? (Now set through constructor)
  int restartFlag {0};
  //PhysicalSystem* physical_system;

  // MPI Communicator for communications among MC walkers
  //MPICommunicator MCAlgorithmCommunicator;

protected :

  double currentTime;
  double lastBackUpTime;

  // it stores the decision of acceptance for each move
  bool acceptMove {false};

  // MC statistics:
  unsigned long int totalMCsteps  {0};
  unsigned long int acceptedMoves {0};
  unsigned long int rejectedMoves {0};

private :


};


// To be implemented.
class Metropolis : public MonteCarloAlgorithm {

public :

  Metropolis(int, const char*) {}
  ~Metropolis() {}

  void run(PhysicalSystem*) {}

};


#endif
