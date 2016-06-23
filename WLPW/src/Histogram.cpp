//#include <cstdlib>
#include <limits>
#include <cmath>
#include <iostream>
#include "Histogram.hpp"

// Constructor
Histogram::Histogram(int restart)
{
  // Markus: Refer to PRE 84, 065702(R) 2011 to set binSize.

  if (restart)
    readHistogramDOSFile("hist_dos_checkpoint.dat");
  else {
    Emax = std::numeric_limits<double>::max();
    //Emin = std::numeric_limits<double>::lowest();    // C++11
    Emin = -std::numeric_limits<double>::max();
    readWangLandauInputFile("wl.input");

    //dim                    = 1;
    //flatnessCriterion      = 0.6;         // change to 0.8 for more accuate results
    //modFactor              = 1.0;         
    //modFactorFinal         = 0.125;       // should be ~ 10E-6
    //modFactorReducer       = 2.0;         // don't change
    //Emin                   = -333.775;      // need to change
    //Emax                   = -333.695;      // need to change
    //binSize                = 0.001;         // need to change
    numBins                = ceil((Emax - Emin) / binSize);
    //histogramCheckInterval = 500;
    //histogramCheckInterval = numBins * 10;

    hist.assign(numBins, 0);
    dos.assign(numBins, 0.0);
    visited.assign(numBins, 0);

    totalMCsteps           = 0;
    acceptedMoves          = 0;
    rejectedMoves          = 0;
    iterations             = 1;
    numBelowRange          = 0;
    numAboveRange          = 0;
  }

  idx           = -1;
  histogramFlat = false;

  printf("Histogram class is created.\n");  
}


// Destructor
Histogram::~Histogram()
{
  //delete[] hist;
  //delete[] dos;
  //delete[] visited;
  hist.clear();
  dos.clear();
  visited.clear();

  printf("Histogram class is destroyed.\n");  

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
  if (idx >= 0)
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
  hist[idx] ++;
  visited[idx] = 1;
  //std::cerr << "idx = " << idx << std::endl;
  //std::cerr << "visited[idx] = " << visited[idx] << std::endl;
}

void Histogram::updateHistogram(double energy)
{
  idx = getIndex(energy);
  hist[idx] ++;
  visited[idx] = 1;
  //std::cerr << "idx = " << idx << std::endl;
  //std::cerr << "visited[idx] = " << visited[idx] << std::endl;
}

void Histogram::updateDOS(double energy)
{
  idx = getIndex(energy);
  dos[idx] += modFactor;
  visited[idx] = 1;
  //std::cerr << "idx = " << idx << std::endl;
  //std::cerr << "visited[idx] = " << visited[idx] << std::endl;
}

bool Histogram::checkEnergyInRange(double energy)
{
  if (energy < Emin) {
    numBelowRange++; 
    return false;
  }
  else if (energy > Emax) {
    numAboveRange++;
    return false;
  }
  else
    return true;
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
  fprintf(histdos_file, "numBelowRange  %ld \n", numBelowRange);
  fprintf(histdos_file, "numAboveRange  %ld \n", numAboveRange);
  
  fprintf(histdos_file, "\n");

  // Write out histogram and DOS
  for (int i=0; i<numBins; i++) {
    fprintf(histdos_file, "%8d %5d %lu %20.8f\n", i, visited[i], hist[i], dos[i]);
    //std::cerr << "visited[" << i << "] = " << visited[i] << std::endl;
  }

  fprintf(histdos_file, "\n");
  fclose(histdos_file);
}


bool Histogram::checkIntegrity()
{

}


// Private member functions
int Histogram::getIndex(double energy)
{
  return floor((energy - Emin) / double(binSize));
}


void Histogram::readHistogramDOSFile(char fileName[])
{

  FILE *histdos_file = fopen(fileName, "r");
  if (histdos_file == NULL)
    std::cerr << "Error: cannot open hist_dos restart file "  << fileName << std::endl;
 
  if (fscanf(histdos_file, "%*s %d", &dim) != 1)
    std::cerr << "Cannot read dim \n";    

  if (fscanf(histdos_file, "%*s %lf", &flatnessCriterion) != 1)
    std::cerr << "Cannot read flatnessCriterion \n";

  if (fscanf(histdos_file, "%*s %lf", &modFactor) != 1)
    std::cerr << "Cannot read modFactor \n";

  if (fscanf(histdos_file, "%*s %lf", &modFactorFinal) != 1)
    std::cerr << "Cannot read modFactorFinal \n";

  if (fscanf(histdos_file, "%*s %lf", &modFactorReducer) != 1)
    std::cerr << "Cannot read modFactorReducer \n";

  if (fscanf(histdos_file, "%*s %lu", &histogramCheckInterval) != 1)
    std::cerr << "Cannot read histogramCheckInterval \n";

  if (fscanf(histdos_file, "%*s %lf", &Emin) != 1)
    std::cerr << "Cannot read Emin \n";

  if (fscanf(histdos_file, "%*s %lf", &Emax) != 1)
    std::cerr << "Cannot read Emax \n";

  if (fscanf(histdos_file, "%*s %lf", &binSize) != 1)
    std::cerr << "Cannot read binSize \n";

  if (fscanf(histdos_file, "%*s %ld", &numBins) != 1)
    std::cerr << "Cannot read numBins \n";

  if (fscanf(histdos_file, "%*s %lu", &totalMCsteps) != 1)
    std::cerr << "Cannot read totalMCsteps \n";

  if (fscanf(histdos_file, "%*s %lu", &acceptedMoves) != 1)
    std::cerr << "Cannot read acceptedMoves \n";

  if (fscanf(histdos_file, "%*s %lu", &rejectedMoves) != 1)
    std::cerr << "Cannot read rejectedMoves \n";

  if (fscanf(histdos_file, "%*s %d", &iterations) != 1)
    std::cerr << "Cannot read iterations \n";

  if (fscanf(histdos_file, "%*s %lu", &numBelowRange) != 1)
    std::cerr << "Cannot read numBelowRange \n";

  if (fscanf(histdos_file, "%*s %lu", &numAboveRange) != 1)
    std::cerr << "Cannot read numAboveRange \n";

  // Open up arrays needed to store the mask, histogram and DOS
  //visited = new int[numBins];
  //hist    = new unsigned long int[numBins];
  //dos     = new double[numBins];
  hist.assign(numBins, 0);
  dos.assign(numBins, 1.0);
  visited.assign(numBins, 0);

  // Continue reading these from file
  int dummy = 0;
  for (int i = 0; i < numBins; i++) {
    if (fscanf(histdos_file, "%d %d %lu %lf", &dummy, &visited[i], &hist[i], &dos[i]) != 4) 
      std::cerr << "Cannot read histogram and DOS.\n";
    if (dummy != i)
      std::cerr << "Problem reading histogram and DOS. Check!\n";
  }

  fclose(histdos_file);

/*
  // YingWai's check, should be removed when things work fine
  std::cerr <<  "YingWai's check for I/O. dim = " << dim << std::endl;
  std::cerr <<  "YingWai's check for I/O. flatnessCriterion = " << flatnessCriterion << std::endl;
  std::cerr <<  "YingWai's check for I/O. modFactor = " << modFactor << std::endl;
  std::cerr <<  "YingWai's check for I/O. modFactorFinal = " << modFactorFinal << std::endl;
  std::cerr <<  "YingWai's check for I/O. modFactorReducer = " << modFactorReducer << std::endl;
  std::cerr <<  "YingWai's check for I/O. histogramCheckInterval = " << histogramCheckInterval << std::endl;
  std::cerr <<  "YingWai's check for I/O. Emin = " << Emin << std::endl;
  std::cerr <<  "YingWai's check for I/O. Emax = " << Emax << std::endl;
  std::cerr <<  "YingWai's check for I/O. binSize = " << binSize << std::endl;
  std::cerr <<  "YingWai's check for I/O. numBins = " << numBins << std::endl;
  std::cerr <<  "\n";

  std::cerr << "YingWai's check for I/O. totalMCsteps = " << totalMCsteps << std::endl;
  std::cerr << "YingWai's check for I/O. acceptedMoves = " << acceptedMoves << std::endl;
  std::cerr << "YingWai's check for I/O. rejectedMoves = " << rejectedMoves << std::endl;
  std::cerr << "YingWai's check for I/O. iterations = " << iterations << std::endl;
  std::cerr << "\n";

  // Write out histogram and DOS
  for (int i = 0; i < numBins; i++) {
    std::cerr << "YingWai's check for I/O. " 
              << i       << " " << visited[i] << " " 
              << hist[i] << " " << dos[i]     << std::endl;
  }
*/

}


void Histogram::readWangLandauInputFile(char fileName[])
{
 
  FILE *WLinput_file = fopen(fileName, "r");
  if (WLinput_file == NULL)
    std::cerr << "Error: cannot open Wang-Landau input file "  << fileName << std::endl;
 
  if (fscanf(WLinput_file, "%*s %d", &dim) != 1)
    std::cerr << "Cannot read dim \n";

  if (fscanf(WLinput_file, "%*s %lf", &flatnessCriterion) != 1)
    std::cerr << "Cannot read flatnessCriterion \n";

  if (fscanf(WLinput_file, "%*s %lf", &modFactor) != 1)
    std::cerr << "Cannot read modFactor \n";

  if (fscanf(WLinput_file, "%*s %lf", &modFactorFinal) != 1)
    std::cerr << "Cannot read modFactorFinal \n";

  if (fscanf(WLinput_file, "%*s %lf", &modFactorReducer) != 1)
    std::cerr << "Cannot read modFactorReducer \n";

  if (fscanf(WLinput_file, "%*s %lu", &histogramCheckInterval) != 1)
    std::cerr << "Cannot read histogramCheckInterval \n";

  if (fscanf(WLinput_file, "%*s %lf", &Emin) != 1)
    std::cerr << "Cannot read Emin \n";

  if (fscanf(WLinput_file, "%*s %lf", &Emax) != 1)
    std::cerr << "Cannot read Emax \n";

  if (fscanf(WLinput_file, "%*s %lf", &binSize) != 1)
    std::cerr << "Cannot read binSize \n";

  fclose(WLinput_file);

}


