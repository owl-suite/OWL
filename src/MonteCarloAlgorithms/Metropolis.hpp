#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP

#include <fstream>
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

  FILE* timeSeriesFile;

  unsigned long int numberOfThermalizationSteps;
  unsigned long int numberOfMCSteps;
  unsigned long int numberOfMCUpdatesPerStep {1};
  double temperature;
  double restartTemperature;

  unsigned long int thermalizationStepsPerformed {0};
  unsigned long int MCStepsPerformed             {0};

  void readMCInputFile(const char* fileName);  // TODO: this should move to MCAlgorithms base class (Histogram class has the same function)
  void readCheckPointFile(const char* fileName);
  void accumulateObservables();                // TODO: this should move to MCAlgorithms base class
  void calculateAveragesAndVariances();

  void writeMCFile(unsigned long int MCSteps);
  void writeStatistics(OutputMode output_mode, const char* = NULL);    // TODO: this should move to MCAlgorithms base class
  void writeCheckPointFiles(OutputMode output_mode);                   // TODO: this should move to MCAlgorithms base class

};

#endif