//#include <cstdlib>
#include <limits>
#include "Histogram.hpp"

//Constructor
Histogram::Histogram()
{
  Dim = 1;
  numBins = 10;
  Emax = std::numeric_limits<double>::max();
  Emin = -std::numeric_limits<double>::infinity();
  //Emin = std::numeric_limits<double>::lowest();    // C++11


  hist = new unsigned long int[numBins];
  dos  = new double[numBins];

  for(int i=0; i<numBins; i++) {
    hist[i] = 0;
    dos[i] = 1.0;
  }

  printf("Histogram class is created.\n");  
}

//Destructor
Histogram::~Histogram()
{
  delete[] hist;
  delete[] dos;
  printf("Histogram class is destroyed.\n");  

}

//Member functions
long int Histogram::getNumberOfBins()
{
  return numBins;
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

void Histogram::setEnergyRange(double E1, double E2)
{
  Emin = E1;
  Emax = E2;
}

