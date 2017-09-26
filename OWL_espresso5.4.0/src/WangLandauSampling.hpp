#ifndef WANG_LANDAU_SAMPLING_HPP
#define WANG_LANDAU_SAMPLING_HPP

#include "MCAlgorithms.hpp"
#include "Histogram.hpp"

class WangLandauSampling : public MonteCarloAlgorithm {

public :

  WangLandauSampling();
  ~WangLandauSampling();

  void run(PhysicalSystem* physical_system);

private :

  Histogram h;

};




#endif
