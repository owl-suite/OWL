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
  logf_final = 0.125;
  Emin       = -340.0;
  Emax       = -330.0;
  binSize    = 0.5;
  numBins    = ceil((Emax - Emin) / binSize);
  numMCSteps = 5;
  //numMCSteps = numBins * 10;

  hist    = new unsigned long int[numBins];
  dos     = new double[numBins];
  visited = new int[numBins];

  for (int i=0; i<numBins; i++)
  {
    hist[i]    = 0;
    dos[i]     = 1.0;
    visited[i] = 0;
  }

  idx = -1;
  totalMCsteps  = 0;
  acceptedMoves = 0;
  rejectedMoves = 0;
  iterations    = 0;

  printf("Histogram class is created.\n");  
}


// Destructor
Histogram::~Histogram()
{
  delete[] hist;
  delete[] dos;
  delete[] visited;
  printf("Histogram class is destroyed.\n");  

}


// Private member functions
int Histogram::getIndex(double energy)
{
  return floor((energy - Emin) / double(binSize));
}


// Public member functions
double Histogram::getBinSize()
{
  return binSize;
}

int Histogram::getNumberOfBins()
{
  return numBins;
}

double Histogram::getDOS(double energy)
{
  idx = getIndex(energy);
  return dos[idx];
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
  // need to resize hist and dos accordingly
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

void Histogram::updateHistogramDOS(double energy)
{
  idx = getIndex(energy);
  dos[idx] += logf;
  hist[idx]++;
  visited[idx] = 1;
}

void Histogram::updateHistogram(double energy)
{
   idx = getIndex(energy);
   hist[idx]++;
   visited[idx] = 1;
}

void Histogram::updateDOS(double energy)
{
  idx = getIndex(energy);
  dos[idx] += logf;
  visited[idx] = 1;
}

bool Histogram::checkHistogramFlatness()
{
  int numVisitedBins = 0;
  unsigned long int sumEntries = 0;
  bool allEntriesPassTest = true;

  for (int i=0; i<numBins; i++) {
    if (visited[i] == 1) {
      sumEntries += hist[i];
      numVisitedBins++;
    }
  }
  double flatnessCriterion = p * double(sumEntries) / double(numVisitedBins);

  for (int i=0; i<numBins; i++) {
    if (visited[i] == 1)
      if (hist[i] < flatnessCriterion) {
        allEntriesPassTest = false;
        break;
      }
  }
  
  return allEntriesPassTest;

}


bool Histogram::checkIntegrity()
{

}
