#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP

#include "MCAlgorithms.hpp"


class Metropolis : public MonteCarloAlgorithm {

public :

  Metropolis(PhysicalSystem* ps, const char* inputFile);
  ~Metropolis();

  void run()                  override;

private :

  PhysicalSystem* physical_system;
  ObservableType* averagedObservables;
  ObservableType* variances;

  FILE* MCOutputFile;

  unsigned long int numberOfThermalizationSteps;
  unsigned long int numberOfMCSteps;
  unsigned long int numberOfMCUpdatesPerStep {1};
  double temperature;

  void readMCInputFile(const char* fileName);  // TODO: this should move to MCAlgorithms base class (Histogram class has the same function)
  void accumulateObservables();                // TODO: this should move to MCAlgorithms base class
  void calculateAveragesAndVariances();
  void writeMCFile(unsigned long int MCSteps);
  void writeResultsFile(const char* = NULL);   // TODO: this should move to MCAlgorithms base class

};

#endif