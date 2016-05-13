//#include <cstdlib>
//#include <limits>
#include <cmath>
#include "Histogram.hpp"

// Constructor
Histogram::Histogram()
{
  // Markus: Refer to PRE 84, 065702(R) 2011 to set binSize.

  //Emax = std::numeric_limits<double>::max();
  //Emin = -std::numeric_limits<double>::infinity();
  //Emin = std::numeric_limits<double>::lowest();    // C++11
  dim                    = 1;
  flatnessCriterion      = 0.6;         // change to 0.8 for more accuate results
  modFactor              = 1.0;         
  modFactorFinal         = 0.125;       // should be ~ 10E-6
  modFactorReducer       = 2.0;         // don't change
  histogramCheckInterval = 5;
  //histogramCheckInterval = numBins * 10;
  Emin                   = -340.0;      // need to change
  Emax                   = -330.0;      // need to change
  binSize                = 0.5;         // need to change
  numBins                = ceil((Emax - Emin) / binSize);

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
  dos[idx] += modFactor;
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
  dos[idx] += modFactor;
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
  double flatnessReference = flatnessCriterion * double(sumEntries) / double(numVisitedBins);

  for (int i=0; i<numBins; i++) {
    if (visited[i] == 1) {
      if (hist[i] < flatnessReference) {
        allEntriesPassTest = false;
        break;
      }
    }
  }
  
  return allEntriesPassTest;

}


bool Histogram::checkIntegrity()
{

}
