#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP


#include <cstdio>
#include <vector>
#include "Main/Globals.hpp"

// TO DO: make it a template class to allow for int / double histogram
class Histogram {

public:

  // These should be moved to the WL-MC class
  // ... or should it be a derived WL histogram class?   (July 24, 2017)
  double       flatnessCriterion;               // flatness criterion  
  double       modFactor;                       // natural log of modification factor, f
  double       modFactorFinal;                  // predefined log(f) to terminate simulation
  double       modFactorReducer;                // a factor to reduce log(f)
  unsigned int histogramCheckInterval;          // number of MC steps between every histogram flatness check
  int          histogramRefreshInterval;        // refresh histogram every certain number of histogram checks

  // MUCA:
  double       KullbackLeiblerDivergence;            // Kullback-Leibler divergence
  double       KullbackLeiblerDivergenceThreshold;   // predefined Kullback-Leibler divergence to terminate simulation
  unsigned int numberOfUpdatesPerIteration;          // number of MC steps in an iteration
  double       numberOfUpdatesMultiplier;            // a multiplier to increase the number of updates in the next iteration
  unsigned int numberOfThermalizationSteps;          // number of MC steps for thermalization before each iteration

  // MC sampling statistics:  (TO DO: these are repeated in MCAlgorithm base class. Revision needed.)
  unsigned long int totalMCsteps;
  unsigned long int acceptedMoves;
  unsigned long int rejectedMoves;
  int  iterations;

  bool histogramFlat;
  int  numHistogramNotImproved;
  int  numHistogramRefreshed;

  // Constructor
  Histogram(int = -1, const char* = NULL, const char* = NULL);
  //Histogram(int = -1, const char* = NULL);

  // Destructor
  ~Histogram();
  
  // Public member functions:
  ObservableType getBinSize();
  unsigned int   getNumberOfBins();
  double         getDOS(ObservableType energy);

  void setEnergyRange (ObservableType E1, ObservableType E2);
  void setBinSize (ObservableType dE);
  void setNumberOfBins (long int n);
  void resetHistogram();
  void refreshHistogram();
  void resetDOS();
  void updateHistogramDOS(ObservableType energy);
  void updateHistogram(ObservableType energy);
  void updateDOS(ObservableType energy);
  void updateDOSwithHistogram();
  void updateDOSwithRemainder();

  void writeHistogramDOSFile(const char* fileName, int iteration = -1, int walkerID = 0);
  void writeNormDOSFile(const char* fileName, int walkerID = 0);

  bool checkEnergyInRange(ObservableType energy);
  bool checkHistogramFlatness();             // for WL
  bool checkKullbackLeiblerDivergence();     // for MUCA
  bool checkIntegrity();                     // check if histogram or DOS have correct bin size,
                                             // number of bins, etc. with respect to the energy range


private:
 
  int dim;                                   // dimension of the histogram

  //ObservableType Emin;                     // energy range for WL sampling (should they be here?)
  double Emin;
  //ObservableType Emax;
  double Emax;
  //ObservableType binSize;                  // energy bin size
  double binSize;
  unsigned int numBins;                      // total number of bins
  unsigned int numVisitedBins;               // number of visited bins
  unsigned int numBinsFailingCriterion;

  unsigned long int numBelowRange;           // count the number of configurations that falls below Emin
  unsigned long int numAboveRange;           // count the number of configurations that falls above Emax

  std::vector<unsigned long int> hist;       // an array to store the histogram
  std::vector<double> dos;                   // an array to store the density of states
  std::vector<int> visited;                  // an array to mark if a bin is visited
  int idx;                                   // index of a bin in the histogram and DOS

  // MUCA only:
  std::vector<double> probDistribution;      // an array to store the probablity distribution constructed from a histogram 

  // REWL only:
  int numberOfWindows;                       // number of energy sub-windows 
  int numberOfWalkersPerWindow;              // number of MC walkers having the same energy sub-windows 
  double overlap;                            // overlapping factor between consecutive energy windows 
  int walkerID;
  int myWindow;

  // Private member functions:
  int getIndex(ObservableType energy);       // Calculate the bin index from an energy
  void readHistogramDOSFile(const char* fileName);
  void readMCInputFile(const char* fileName);
  void calculateProbabilityDistribution();   // Normalized histogram
  unsigned int getNumberOfVisitedBins();
  void shiftDOS();

};

#endif

