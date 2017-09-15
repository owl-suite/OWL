#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP


#include <cstdio>
#include <vector>

//typedef int ObservableType;
typedef double ObservableType;


// TO DO: make it a template class to allow for int / double histogram
class Histogram {

public:

  // These should be moved to the WL-MC class (when they are implemented...)
  // ... or should it be a derived WL histogram class?   (July 24, 2017)
  double       flatnessCriterion;               // flatness criterion  
  double       modFactor;                       // natural log of modification factor, f
  double       modFactorFinal;                  // predefined log(f) to terminate simulation
  double       modFactorReducer;                // a factor to reduce log(f)
  unsigned int histogramCheckInterval;          // number of MC steps between every histogram flatness check
  int          histogramRefreshInterval;        // refresh histogram every certain number of histogram checks

  // MUCA statistics:
  double       KullbackLeiblerDivergence;
  double       KullbackLeiblerDivergenceThreshold;

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
  ObservableType getBinSize();
  int            getNumberOfBins();
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

  void writeHistogramDOSFile(const char* fileName);
  void writeNormDOSFile(const char* fileName);

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
  unsigned int numBinsFailingCriterion;

  unsigned long int numBelowRange;           // count the number of configurations that falls below Emin
  unsigned long int numAboveRange;           // count the number of configurations that falls above Emax

  std::vector<unsigned long int> hist;       // an array to store the histogram
  std::vector<double> dos;                   // an array to store the density of states
  std::vector<int> visited;                  // an array to mark if a bin is visited
  int idx;                                   // index of a bin in the histogram and DOS

  std::vector<double> probDistribution;      // an array to store the probablity distribution constructed from a histogram  (MUCA only)

  // Private member functions:
  int getIndex(ObservableType energy);       // Calculate the bin index from an energy
  void readHistogramDOSFile(const char* fileName);
  void readWangLandauInputFile(const char* fileName);

};

#endif

