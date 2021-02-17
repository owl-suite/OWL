#ifndef GLOBAL_UPDATE_WANG_LANDAU_SAMPLING_HPP
#define GLOBAL_UPDATE_WANG_LANDAU_SAMPLING_HPP

#include "MCAlgorithms.hpp"
#include "Histogram.hpp"

class GlobalUpdateWangLandauSampling : public MonteCarloAlgorithm {

public :

  GlobalUpdateWangLandauSampling(PhysicalSystem* ps);
  ~GlobalUpdateWangLandauSampling();

  void run()                  override;

private :

  PhysicalSystem* physical_system;
  Histogram h;
  
};

#endif
