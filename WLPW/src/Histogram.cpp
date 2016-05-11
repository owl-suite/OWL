//#include <cstdlib>
//#include <limits>
#include <cmath>
#include "Histogram.hpp"

// Constructor
Histogram::Histogram()
{
  //Emax = std::numeric_limits<double>::max();
  //Emin = -std::numeric_limits<double>::infinity();
  //Emin = std::numeric_limits<double>::lowest();    // C++11
  dim        = 1;
  p          = 0.6;
  logf       = 1.0;
  logf_final = 1E-6;
  Emin       = -400.0;
  Emax       = -300.0;
  binSize    = 1.0;
  numBins    = ceil((Emax - Emin) / binSize);
  numMCSteps = numBins * 10;

  hist = new unsigned long int[numBins];
  dos  = new double[numBins];

  for (int i=0; i<numBins; i++)
  {
    hist[i] = 0;
    dos[i] = 1.0;
  }

  totalMCsteps = 0;
  nIteration = 0;

  printf("Histogram class is created.\n");  
}


// Destructor
Histogram::~Histogram()
{
  delete[] hist;
  delete[] dos;
  printf("Histogram class is destroyed.\n");  

}


//Member functions
double Histogram::getBinSize()
{
  return binSize;
}

int Histogram::getNumberOfBins()
{
  return numBins;
}

void Histogram::setEnergyRange(double E1, double E2)
{
  Emin = E1;
  Emax = E2;
}

void Histogram::setBinSize(double dE)
{
  binSize = dE;
}

void Histogram::setNumberOfBins(long int n)
{
  numBins = n;
}

void Histogram::resetHistogram()
{
  for (int i=0; i<numBins; i++)
    hist[i] = 0;
}

void Histogram::resetDOS()
{
  for (int i=0; i<numBins; i++)
    dos[i] = 1.0;
}

bool Histogram::checkIntegrity()
{

}
