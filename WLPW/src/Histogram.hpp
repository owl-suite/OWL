#include <cstdio>

class Histogram {

private:
 
  // WL sampling parameters:
  int dim;                       // dimension of the histogram

  double p;                      // flatness criterion  
  double logf;                   // natural log of modification factor
  double logf_final;             // predefined log(f) to terminate simulation

  double Emin;                   // energy range for WL sampling (should they be here?)
  double Emax;
  double binSize;                // energy bin size
  int numBins;                   // total number of bins
  int numMCSteps;                // number of MC steps between every histogram flatness check

  unsigned long int *hist;       // an array to store the histogram
  double *dos;                   // an array to store the density of states

  // WL sampling statistics:
  unsigned long int totalMCsteps;
  int nIteration;

public:
  //Constructor
  Histogram();
  //Destructor
  ~Histogram();
  
  //Member functions:
  double getBinSize();
  int getNumberOfBins();

  void setEnergyRange (double, double);
  void setBinSize (double);
  void setNumberOfBins (long int);
  void resetHistogram();
  void resetDOS();

  bool checkIntegrity();          // check if histogram or DOS have correct bin size,
                                  // number of bins, etc. with respect to the energy range

  
};
