#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP


#include <cstdio>
#include <vector>

class Histogram {

public:

  // These should be moved to the WL-MC class (when they are implemented...)
  double       flatnessCriterion;               // flatness criterion  
  double       modFactor;                       // natural log of modification factor, f
  double       modFactorFinal;                  // predefined log(f) to terminate simulation
  double       modFactorReducer;                // a factor to reduce log(f)
  unsigned int histogramCheckInterval;          // number of MC steps between every histogram flatness check
  int          histogramRefreshInterval;        // refresh histogram every certain number of histogram checks

  // WL sampling statistics:
  unsigned long int totalMCsteps;
  unsigned long int acceptedMoves;
  unsigned long int rejectedMoves;
  int  iterations;

  bool histogramFlat;
  int  numHistogramNotImproved;
  int  numHistogramRefreshed;

  // Constructor
  Histogram(int = -1, const char* = NULL);

  // Destructor
  ~Histogram();
  
  // Public member functions:
  double getBinSize();
  int    getNumberOfBins();
  double getDOS(double);

  void setEnergyRange (double, double);
  void setBinSize (double);
  void setNumberOfBins (long int);
  void resetHistogram();
  void refreshHistogram();
  void resetDOS();
  void updateHistogramDOS(double);
  void updateHistogram(double);
  void updateDOS(double);

  void writeHistogramDOSFile(const char*);
  void writeNormDOSFile(const char*);

  bool checkEnergyInRange(double energy);
  bool checkHistogramFlatness();
  bool checkIntegrity();          // check if histogram or DOS have correct bin size,
                                  // number of bins, etc. with respect to the energy range


private:
 
  int dim;                                   // dimension of the histogram

  double Emin;                               // energy range for WL sampling (should they be here?)
  double Emax;
  double binSize;                            // energy bin size
  unsigned int numBins;                      // total number of bins
  unsigned int numBinsFailingCriterion;

  unsigned long int numBelowRange;           // count the number of configurations that falls below Emin
  unsigned long int numAboveRange;           // count the number of configurations that falls above Emax

  std::vector<unsigned long int> hist;       // an array to store the histogram
  std::vector<double> dos;                   // an array to store the density of states
  std::vector<int> visited;                  // an array to mark if a bin is visited
  int idx;                                   // index of a bin in the histogram and DOS

  // Private member functions:
  int getIndex(double);                      // Calculate the bin index from an energy
  void readHistogramDOSFile(char const[]);
  void readWangLandauInputFile(char const[]);

};

#endif

