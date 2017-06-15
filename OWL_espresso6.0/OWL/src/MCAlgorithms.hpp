#ifndef MC_ALGORITHMS_HPP
#define MC_ALGORITHMS_HPP

#include <mpi.h>
#include "PhysicalSystemBase.hpp"
#include "Histogram.hpp"
//# include "Communications.hpp"


// Base class for all Monte Carlo algorithms
class MonteCarloAlgorithm {

public :

  // Constructor
  MonteCarloAlgorithm();

  // Destructor
  virtual ~MonteCarloAlgorithm() {}

  virtual void run(PhysicalSystem*) = 0;

  // Should they be set directly by I/O? (Now set through constructor)
  int restartFlag;     
  //PhysicalSystem* physical_system;

protected :

  double currentTime;
  double lastBackUpTime;

  // it stores the decision of WL acceptance for each move
  bool acceptMove {false};

  // MC statistics:
  unsigned long int totalMCsteps;
  unsigned long int acceptedMoves;
  unsigned long int rejectedMoves;

private :


};



class WangLandauSampling : public MonteCarloAlgorithm {

public :
  
  WangLandauSampling(int, const char*);
  ~WangLandauSampling();

  void run(PhysicalSystem*);

private :

  Histogram h;

};


// To be implemented.
class Metropolis : public MonteCarloAlgorithm {

public :

  Metropolis(int, const char*) {}
  ~Metropolis() {}

  void run(PhysicalSystem*) {}

};


#endif
