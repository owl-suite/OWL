#ifndef MULTICANONICAL_SAMPLING_HPP
#define MULTICANONICAL_SAMPLING_HPP

#include "MCAlgorithms.hpp"
#include "Histogram.hpp"

class MulticanonicalSampling : public MonteCarloAlgorithm {

public :

  MulticanonicalSampling(PhysicalSystem* ps);
  ~MulticanonicalSampling();

  void run();

private :

  PhysicalSystem* physical_system;
  Histogram h;

};




#endif
