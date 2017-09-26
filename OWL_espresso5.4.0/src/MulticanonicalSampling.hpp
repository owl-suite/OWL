#ifndef MULTICANONICAL_SAMPLING_HPP
#define MULTICANONICAL_SAMPLING_HPP

#include "MCAlgorithms.hpp"
#include "Histogram.hpp"

class MulticanonicalSampling : public MonteCarloAlgorithm {

public :

  MulticanonicalSampling();
  ~MulticanonicalSampling();

  void run(PhysicalSystem* physical_system);

private :

  Histogram h;

};




#endif
