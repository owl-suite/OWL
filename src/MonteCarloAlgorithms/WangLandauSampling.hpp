#ifndef WANG_LANDAU_SAMPLING_HPP
#define WANG_LANDAU_SAMPLING_HPP

#include "MCAlgorithms.hpp"
#include "Histogram.hpp"

class WangLandauSampling : public MonteCarloAlgorithm {

public :

  WangLandauSampling(PhysicalSystem* ps);
  ~WangLandauSampling();

  void run();

private :

  PhysicalSystem* physical_system;
  Histogram h;
  
};

#endif
