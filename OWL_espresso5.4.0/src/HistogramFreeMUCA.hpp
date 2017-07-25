#ifndef HISTOGRAM_FREE_MUCA_HPP
#define HISTOGRAM_FREE_MUCA_HPP

#include <vector>
#include "MCAlgorithms.hpp"
#include "Histogram.hpp"


// TO DO: should be changed into a template to allow for int / double DataSet  (July 16, 2017)
class DiscreteHistogramFreeMUCA : public MonteCarloAlgorithm {

public :

  DiscreteHistogramFreeMUCA(int restart, const char* inputFile);
  ~DiscreteHistogramFreeMUCA();

  void run(PhysicalSystem*);

private :

  Histogram h;                  // ironically a histogram is still needed for the discrete case
  int numberOfDataPoints;       // number of data points in each data set
  std::vector<int> DataSet;     // The list of energy (data set) in each iteration

  void resetDataSet();

};




#endif
