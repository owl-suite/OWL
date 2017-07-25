#ifndef MULTICANONICAL_SAMPLING_HPP
#define MULTICANONICAL_SAMPLING_HPP

#include "MCAlgorithms.hpp"
#include "Histogram.hpp"

class MulticanonicalSampling : public MonteCarloAlgorithm {

public :

  MulticanonicalSampling(int restart, const char* inputFile);
  ~MulticanonicalSampling();

  void run(PhysicalSystem*);

private :

  Histogram h;

};




#endif
