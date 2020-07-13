// MUCA implementation following the recipe in this paper:
// Ref: J. Gross, J. Zierenberg, M. Weigel, and W. Janke. Comp. Phys. Comm. 224, 387–395 (2018). 

#ifndef MULTICANONICAL_SAMPLING_HPP
#define MULTICANONICAL_SAMPLING_HPP

#include "MCAlgorithms.hpp"
#include "Histogram.hpp"

class MulticanonicalSampling : public MonteCarloAlgorithm {

public :

  MulticanonicalSampling(PhysicalSystem* ps);
  ~MulticanonicalSampling();

  void run()                  override;

private :

  PhysicalSystem* physical_system;
  Histogram h;

};




#endif
