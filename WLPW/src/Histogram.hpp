#include <cstdio>

class Histogram {

private:
  int Dim;                           // dimension of the histogram
  long int numBins;                  // total number of bins
  unsigned long int *hist;           // an array to store the histogram
  double *dos;                       // an array to store the density of states
  double Emin;                       // energy range for WL sampling (should they be here?)
  double Emax;


public:
  //Constructor
  Histogram();
  //Destructor
  ~Histogram();
  
  //Member functions:
  long int getNumberOfBins();
  void setNumberOfBins(long int);
  void resetHistogram();
  void resetDOS();
  void setEnergyRange(double, double);

};
