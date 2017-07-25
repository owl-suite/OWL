#ifndef WANG_LANDAU_SAMPLING_HPP
#define WANG_LANDAU_SAMPLING_HPP

#include "MCAlgorithms.hpp"
#include "Histogram.hpp"

class WangLandauSampling : public MonteCarloAlgorithm {

public :

  WangLandauSampling(int restart, const char* inputFile);
  ~WangLandauSampling();

  void run(PhysicalSystem*);

private :

  Histogram h;

};




#endif
