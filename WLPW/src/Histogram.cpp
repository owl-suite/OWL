//#include <cstdlib>
//#include <limits>
#include <cmath>
#include "Histogram.hpp"

// Constructor #1 : for fresh run
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
  Emin                   = -334.0;      // need to change
  Emax                   = -333.0;      // need to change
  binSize                = 0.01;        // need to change
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
  histogramFlat = false;

  printf("Histogram class is created.\n");  
}


// Constructor #2 : for restarted run
Histogram::Histogram(char fileName[])
{
  readHistogramDOSFile(fileName);

  idx = -1;
  histogramFlat = false;

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

void Histogram::writeHistogramDOSFile(char fileName[])
{
  FILE *histdos_file;
  histdos_file = fopen(fileName, "w");

  // Write out histogram info
  fprintf(histdos_file, "dim  %d \n", dim);
  fprintf(histdos_file, "flatnessCriterion  %5.3f \n", flatnessCriterion);
  fprintf(histdos_file, "modFactor  %5.3f \n", modFactor);
  fprintf(histdos_file, "modFactorFinal  %12.8e \n", modFactorFinal);
  fprintf(histdos_file, "modFactorReducer  %5.3f \n", modFactorReducer);
  fprintf(histdos_file, "histogramCheckInterval  %ld \n", histogramCheckInterval);
  fprintf(histdos_file, "Emin  %15.8e \n", Emin);
  fprintf(histdos_file, "Emax  %15.8e \n", Emax);
  fprintf(histdos_file, "binSize  %15.8e \n", binSize);
  fprintf(histdos_file, "numBins  %ld \n", numBins);
  fprintf(histdos_file, "\n");

  // Write out WL statistics
  fprintf(histdos_file, "totalMCsteps  %ld \n", totalMCsteps);
  fprintf(histdos_file, "acceptedMoves  %ld \n", acceptedMoves);
  fprintf(histdos_file, "rejectedMoves  %ld \n", rejectedMoves);
  fprintf(histdos_file, "iterations  %d \n", iterations);
  fprintf(histdos_file, "\n");

  // Write out histogram and DOS
  for (int i=0; i<numBins; i++)
    fprintf(histdos_file, "%8d %5d %20d %20.8f\n", i, visited[i], hist[i], dos[i]);

  fprintf(histdos_file, "\n");
  fclose(histdos_file);
}

void Histogram::readHistogramDOSFile(char fileName[])
{

  // read these from file:
  dim                    = 1;
  flatnessCriterion      = 0.6;         // change to 0.8 for more accuate results
  modFactor              = 1.0;         
  modFactorFinal         = 0.125;       // should be ~ 10E-6
  modFactorReducer       = 2.0;         // don't change
  histogramCheckInterval = 5;
  Emin                   = -340.0;      // need to change
  Emax                   = -330.0;      // need to change
  binSize                = 0.5;         // need to change
  numBins                = ceil((Emax - Emin) / binSize);
  totalMCsteps  = 0;
  acceptedMoves = 0;
  rejectedMoves = 0;
  iterations    = 0;

  // Open up arrays needed to store the mask, histogram and DOS
  hist    = new unsigned long int[numBins];
  dos     = new double[numBins];
  visited = new int[numBins];

  // Continue reading these from file
  for (int i=0; i<numBins; i++)
  {
    hist[i]    = 0;
    dos[i]     = 1.0;
    visited[i] = 0;
  }

}


bool Histogram::checkIntegrity()
{

}
