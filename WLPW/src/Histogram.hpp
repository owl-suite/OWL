#include <cstdio>

class Histogram {

private:
 
  int dim;                       // dimension of the histogram

  double Emin;                   // energy range for WL sampling (should they be here?)
  double Emax;
  double binSize;                // energy bin size
  int numBins;                   // total number of bins

  unsigned long int *hist;       // an array to store the histogram
  double *dos;                   // an array to store the density of states
  int *visited;                  // an array to mark if a bin is visited

  int idx;                       // index of a bin in the histogram and DOS


  // Private member functions:
  int getIndex(double);          // Calculate the bin index from an energy

public:

  double p;                      // flatness criterion  
  double logf;                   // natural log of modification factor
  double logf_final;             // predefined log(f) to terminate simulation
  int numMCSteps;                // number of MC steps between every histogram flatness check

  // WL sampling statistics:
  unsigned long int totalMCsteps;
  unsigned long int acceptedMoves;
  unsigned long int rejectedMoves;
  int iterations;
  bool histogramFlat; 

 
  // Constructor
  Histogram();
  // Destructor
  ~Histogram();
  
  // Public member functions:
  double getBinSize();
  int getNumberOfBins();
  double getDOS(double);

  void setEnergyRange (double, double);
  void setBinSize (double);
  void setNumberOfBins (long int);
  void resetHistogram();
  void resetDOS();
  void updateHistogramDOS(double);
  void updateHistogram(double);
  void updateDOS(double);

  bool checkHistogramFlatness();
  bool checkIntegrity();          // check if histogram or DOS have correct bin size,
                                  // number of bins, etc. with respect to the energy range

  
};
