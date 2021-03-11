#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>             // std::string
#include <sstream>            // std::istringstream
#include "Histogram.hpp"
#include "Main/Communications.hpp"
#include "Utilities/CheckFile.hpp"

//TODO: The argument list of the Histogram constructor is redundant as simInfo is a global struct. Should be removed. (Feb 1, 19)
//TODO: Markus: Refer to PRE 84, 065702(R) 2011 to set binSize.

// Constructor
Histogram::Histogram(int restart, const char* inputFile, const char* checkPointFile)
{

  std::cout << "\nInitializing histogram...\n";

  Emax = std::numeric_limits<ObservableType>::max();
  Emin = std::numeric_limits<ObservableType>::lowest();

  numberOfWindows          = 1;           // default to =1 if not specified in input file
  numberOfWalkersPerWindow = 1;
  overlap  = 1.0;
  walkerID = 0;
  myWindow = 0;

  // Read input file
  if ( file_exists(inputFile) )
    readMCInputFile(inputFile);
  else {
    std::cout << "   Error: No input file for reading histogram's info. Quiting... \n";
    exit(7);
  }

  // Read checkpoint file, or initialize anew
  if (restart && file_exists(checkPointFile) )
    readHistogramDOSFile(checkPointFile);
  else {
    // Calculate quantities based on the variables read in
    double energySubwindowWidth = (Emax - Emin) / (1.0 + double(numberOfWindows - 1)*(1.0 - overlap));

    // Round upward to the closest binSize
    energySubwindowWidth = ceil(energySubwindowWidth / binSize) * binSize;

    walkerID = (GlobalComm.thisMPIrank - (GlobalComm.thisMPIrank % simInfo.numMPIranksPerWalker)) / simInfo.numMPIranksPerWalker;
    myWindow = ( walkerID - (walkerID % numberOfWalkersPerWindow) ) / numberOfWalkersPerWindow;
    Emin     = Emin + double(myWindow) * (1.0 - overlap) * energySubwindowWidth;
    Emax     = Emin + energySubwindowWidth;
    numBins  = unsigned(ceil((Emax - Emin) / binSize)) + 1;

    // YingWai's check
    //printf("YingWai's check: Inside Histogram constructor. world_rank = %3d, myWindow = %3d, walkerID = %3d, Emin = %6.3f, Emax = %6.3f \n", GlobalComm.thisMPIrank, myWindow, walkerID, Emin, Emax);

    hist.assign(numBins, 0);
    dos.assign(numBins, 0.0);
    visited.assign(numBins, 0);
    probDistribution.assign(numBins, 0.0);

    totalMCsteps            = 0;
    acceptedMoves           = 0;
    rejectedMoves           = 0;
    iterations              = 0;
    numBinsFailingCriterion = numBins;
    numBelowRange           = 0;
    numAboveRange           = 0;
    numHistogramNotImproved = 0;
    numHistogramRefreshed   = 0;

    // MUCA:
    KullbackLeiblerDivergence = 0.0;

  }

  idx           = -1;
  histogramFlat = false;

  //printf("Histogram class is created.\n");  
}


// Destructor
Histogram::~Histogram()
{
  hist.clear();
  dos.clear();
  visited.clear();

  //printf("Histogram class is destroyed.\n");  

}


// Public member functions
ObservableType Histogram::getBinSize()
{
  return binSize;
}


unsigned int Histogram::getNumberOfBins()
{
  return numBins;
}


double Histogram::getDOS(ObservableType energy)
{
  idx = getIndex(energy);
  if (idx >= 0)
    return dos[unsigned(idx)];
  else {
    std::cerr << "Problem in File " << __FILE__ << " Line " << __LINE__  << "\n";
    exit(EXIT_FAILURE);
  }
  
}


void Histogram::setEnergyRange(ObservableType E1, ObservableType E2)
{
  Emin = E1;
  Emax = E2;
}


void Histogram::setBinSize(ObservableType dE)
{
  binSize = dE;
}


void Histogram::setNumberOfBins(unsigned int n)
{
  numBins = n;
  // need to resize hist and dos accordingly
}


void Histogram::resetHistogram()
{
  for (unsigned int i=0; i<numBins; i++) {
    hist[i] = 0;
    //std::fill(hist.begin(), hist.end(), 0);      // Intel compiler has problem with it!
    probDistribution[i] = 0.0;
  }
  numHistogramNotImproved = 0;
  numBinsFailingCriterion = numBins;
  numHistogramRefreshed = 0;
}


void Histogram::refreshHistogram()
{
  for (unsigned int i=0; i<numBins; i++) {
    hist[i] = 0;
    probDistribution[i] = 0.0;
  }

  numHistogramNotImproved = 0;
  numBinsFailingCriterion = numBins;
  numHistogramRefreshed++;
}


void Histogram::resetDOS()
{
  for (unsigned int i=0; i<numBins; i++)
    dos[i] = 0.0;
    //std::fill(dos.begin(), dos.end(), 0.0);      // Intel compiler has problem with it!
}


void Histogram::updateHistogramDOS(ObservableType energy)
{
  idx = getIndex(energy);

  if ( idx >= 0 ) {
    unsigned int index = unsigned(idx);
    // If it is the first time a bin is visited:
    //   1. see if it can reference the DOS from neighboring bins
    //   2. reset Histogram and start over
    if ( visited[index] == 0 ) {
      unsigned int refIdx = index;
      if ( index == 0 )
        refIdx = 1;
      else if ( index == (numBins-1) )
        refIdx = index - 1;
      else {
        if ( (visited[index-1] > 0) && (visited[index+1] > 0) )
          refIdx = ( dos[index-1] < dos[index+1] ? (index-1) : (index+1) );
        else if ( visited[index-1] == 0 )
          refIdx = index + 1;
        else if ( visited[index+1] == 0 )
          refIdx = index - 1;
      }

      dos[index] = dos[refIdx];
      visited[index] = 1;
      refreshHistogram();
    }
    else {
      dos[index] += modFactor;
      hist[index]++;
      //visited[index] = 1;
    }
  }
  else { 
    std::cerr << "Error: idx < 0 in updateHistogramDOS!!\n";
    std::cerr << "Aborting...\n";
    exit(10);
  }
  //std::cerr << "energy = " << energy << "\n";
  //std::cerr << "idx = " << idx << "\n";
  //std::cerr << "visited[idx] = " << visited[idx] << "\n";
  //std::cerr << "hist[idx] = " << hist[idx] << "\n";
}


void Histogram::globalUpdateHistogramDOS(ObservableType energy)
{
  idx = getIndex(energy);

  if ( idx >= 0 ) {
    unsigned int index = unsigned(idx);
    // If it is the first time a bin is visited:
    //   1. see if it can reference the DOS from neighboring bins
    //   2. reset Histogram and start over
    if ( visited[index] == 0 ) {
      unsigned int refIdx = index;
      if ( index == 0 )
        refIdx = 1;
      else if ( index == (numBins-1) )
        refIdx = index - 1;
      else {
        if ( (visited[index-1] > 0) && (visited[index+1] > 0) )
          refIdx = ( dos[index-1] < dos[index+1] ? (index-1) : (index+1) );
        else if ( visited[index-1] == 0 )
          refIdx = index + 1;
        else if ( visited[index+1] == 0 )
          refIdx = index - 1;
      }

      dos[index] = dos[refIdx];
      visited[index] = 1;
      refreshHistogram();
    }
    else {
//      dos[index]++;
//      hist[index]++;
//      for (unsigned int i=0; i<numBins; i++)
//        dos[i] -= 1.0 / double(numBins);
      dos[index] += modFactor;
      hist[index]++;
      for (unsigned int i=0; i<numBins; i++)
        dos[i] -= modFactor / double(numBins);
    }
  }
  else { 
    std::cerr << "Error: idx < 0 in globalUpdateHistogramDOS!!\n";
    std::cerr << "Aborting...\n";
    exit(10);
  }
  //std::cerr << "energy = " << energy << "\n";
  //std::cerr << "idx = " << idx << "\n";
  //std::cerr << "visited[idx] = " << visited[idx] << "\n";
  //std::cerr << "hist[idx] = " << hist[idx] << "\n";
}


void Histogram::updateHistogram(ObservableType energy)
{
  idx = getIndex(energy);
  hist[unsigned(idx)]++;
  visited[unsigned(idx)] = 1;
  //std::cerr << "idx = " << idx << "\n";
  //std::cerr << "visited[idx] = " << visited[idx] << "\n";
}

void Histogram::updateDOS(ObservableType energy)
{
  idx = getIndex(energy);
  dos[unsigned(idx)] += modFactor;
  visited[unsigned(idx)] = 1;
  //std::cerr << "idx = " << idx << "\n";
  //std::cerr << "visited[idx] = " << visited[idx] << "\n";
}

void Histogram::updateDOSwithHistogram()
{
  for (unsigned int i=0; i<numBins; i++)
    if ((visited[i] == 1) && (hist[i] != 0))
      dos[i] += log(double(hist[i]));
}

bool Histogram::checkEnergyInRange(ObservableType energy)
{
  bool isWithinRange {false};
  if (energy < Emin) {
    //std::cerr << "Energy below range. Energy = " << energy << "\n";
    numBelowRange++; 
    isWithinRange = false;
  }
  else if (energy > Emax) {
    //std::cerr << "Energy above range. Energy = " << energy << "\n";
    numAboveRange++;
    isWithinRange = false;
  }
  else {
    //std::cerr << "Energy within range. Energy = " << energy << "\n";
    isWithinRange = true;
  }

  return isWithinRange;

}

bool Histogram::checkHistogramFlatness()
{
  int numVisitedBins = 0;
  unsigned long int sumEntries = 0;

  for (unsigned int i=0; i<numBins; i++) {
    if (visited[i] == 1) {
      sumEntries += hist[i];
      numVisitedBins++;
    }
  }
  double flatnessReference = flatnessCriterion * double(sumEntries) / double(numVisitedBins);

  unsigned int numBinsFailingCriterionLastTime = numBinsFailingCriterion;
  numBinsFailingCriterion = 0;
  for (unsigned int i=0; i<numBins; i++) {
    if (visited[i] == 1) {
      if (double(hist[i]) < flatnessReference) {
        numBinsFailingCriterion++;
      } 
    }
  }

  if (numBinsFailingCriterion >= numBinsFailingCriterionLastTime)
    numHistogramNotImproved++;

  if (numBinsFailingCriterion == 0)
    return true;
  else
    return false;

}


// Ref: S. Kullback and R. A. Leibler, Ann. Math. Stat. 22, 79 (1951).
// It measures the similarity of two probability distributions, P(x) and Q(x).
bool Histogram::checkKullbackLeiblerDivergence()
{
  long int numVisitedBins = std::count(visited.begin(), visited.end(), 1);
  //int numVisitedBins = 0;
  //for (unsigned int i=0; i<numBins; i++) {
  //  if (visited[i] == 1)
  //    numVisitedBins++;
  //}
  std::cout << "Number of visited bins = " << numVisitedBins << "\n";

  double flatnessReference = 1.0 / static_cast<double>( numVisitedBins );
  //double flatnessReference = 1.0 / static_cast<double>( std::max(numVisitedBins, 10) );

  KullbackLeiblerDivergence = 0.0;
  for (unsigned int i=0; i<numBins; i++) {
    if (visited[i] == 1) {
      probDistribution[i] = static_cast<double>(hist[i]) / static_cast<double>(numberOfUpdatesPerIteration);
      KullbackLeiblerDivergence += probDistribution[i] * log(probDistribution[i]/flatnessReference);
    }
  }

  std::cout << "KullbackLeiblerDivergence = " << KullbackLeiblerDivergence << "\n";

  if (KullbackLeiblerDivergence <= KullbackLeiblerDivergenceThreshold)
    return true;
  else  
    return false;

}


void Histogram::writeHistogramDOSFile(const char* fileName)
{

  FILE *histdos_file;
  histdos_file = fopen(fileName, "w");

  // Write out histogram info
  fprintf(histdos_file, "dim  %d \n", dim);
  fprintf(histdos_file, "flatnessCriterion  %5.3f \n", flatnessCriterion);
  fprintf(histdos_file, "modFactor  %12.8e \n", modFactor);
  fprintf(histdos_file, "modFactorFinal  %12.8e \n", modFactorFinal);
  fprintf(histdos_file, "modFactorReducer  %5.3f \n", modFactorReducer);
  fprintf(histdos_file, "histogramCheckInterval  %10u \n", histogramCheckInterval);
  fprintf(histdos_file, "histogramRefreshInterval  %10u \n", histogramRefreshInterval);
  fprintf(histdos_file, "Emin  %15.8e \n", Emin);
  fprintf(histdos_file, "Emax  %15.8e \n", Emax);
  fprintf(histdos_file, "binSize  %15.8e \n", binSize);
  fprintf(histdos_file, "numBins  %8u \n", numBins);
  fprintf(histdos_file, "\n");

  // Write out WL statistics
  fprintf(histdos_file, "totalMCsteps  %ld \n", totalMCsteps);
  fprintf(histdos_file, "acceptedMoves  %ld \n", acceptedMoves);
  fprintf(histdos_file, "rejectedMoves  %ld \n", rejectedMoves);
  fprintf(histdos_file, "iterations  %d \n", iterations);
  fprintf(histdos_file, "numBelowRange  %ld \n", numBelowRange);
  fprintf(histdos_file, "numAboveRange  %ld \n", numAboveRange);
  fprintf(histdos_file, "numBinsFailingCriterion  %8u \n", numBinsFailingCriterion);
  fprintf(histdos_file, "numHistogramNotImproved  %d \n", numHistogramNotImproved);
  fprintf(histdos_file, "numHistogramRefreshed %d \n", numHistogramRefreshed);
 
  fprintf(histdos_file, "\n");

  // Write out histogram and DOS
  for (unsigned int i=0; i<numBins; i++) {
    fprintf(histdos_file, "%8d %5d %lu %20.8f\n", i, visited[i], hist[i], dos[i]);
    //std::cerr << "visited[" << i << "] = " << visited[i] << "\n";
  }

  fprintf(histdos_file, "\n");
  fclose(histdos_file);
}


bool Histogram::checkIntegrity()
{
  return true;
}


// Private member functions
int Histogram::getIndex(ObservableType energy)
{
  return int( floor(double(energy - Emin) / double(binSize)) );
}


void Histogram::readHistogramDOSFile(const char* fileName)
{

  if (GlobalComm.thisMPIrank == 0)
    std::cout << "   Reading histogram checkpoint file : " << fileName << "\n";

  FILE *histdos_file = fopen(fileName, "r");
  if (histdos_file == NULL) {
    std::cerr << "     ERROR! Cannot open histogram checkpoint file "  << fileName << "\n";
    exit(7);    // perhaps can start from wl.input instead of quitting?
  }

  // TO DO: need to error-proof if they are not in order / missing...
  if (fscanf(histdos_file, "%*s %d", &dim) != 1)
    std::cerr << "     ERROR! Cannot read dim \n";    

  if (fscanf(histdos_file, "%*s %lf", &flatnessCriterion) != 1)
    std::cerr << "     ERROR! Cannot read flatnessCriterion \n";

  if (fscanf(histdos_file, "%*s %lf", &modFactor) != 1)
    std::cerr << "     ERROR! Cannot read modFactor \n";

  double modFactorFinal_tmp;
  if (fscanf(histdos_file, "%*s %lf", &modFactorFinal_tmp) != 1)
    std::cerr << "     ERROR! Cannot read modFactorFinal \n";
  if (modFactorFinal_tmp < modFactorFinal)
    std::cout << "     WARNING! modFactorFinal read from the checkpoint file is smaller than the input file. Simulation will continue with the new modFactorFinal. \n";

  if (fscanf(histdos_file, "%*s %lf", &modFactorReducer) != 1)
    std::cerr << "     ERROR! Cannot read modFactorReducer \n";

  if (fscanf(histdos_file, "%*s %u", &histogramCheckInterval) != 1)
    std::cerr << "     ERROR! Cannot read histogramCheckInterval \n";

  if (fscanf(histdos_file, "%*s %u", &histogramRefreshInterval) != 1)
    std::cerr << "     ERROR! Cannot read histogramRefreshInterval \n";

  if (fscanf(histdos_file, "%*s %lf", &Emin) != 1)
    std::cerr << "     ERROR! Cannot read Emin \n";

  if (fscanf(histdos_file, "%*s %lf", &Emax) != 1)
    std::cerr << "     ERROR! Cannot read Emax \n";

  if (fscanf(histdos_file, "%*s %lf", &binSize) != 1)
    std::cerr << "     ERROR! Cannot read binSize \n";

  if (fscanf(histdos_file, "%*s %u", &numBins) != 1)
    std::cerr << "     ERROR! Cannot read numBins \n";

  if (fscanf(histdos_file, "%*s %lu", &totalMCsteps) != 1)
    std::cerr << "     ERROR! Cannot read totalMCsteps \n";

  if (fscanf(histdos_file, "%*s %lu", &acceptedMoves) != 1)
    std::cerr << "     ERROR! Cannot read acceptedMoves \n";

  if (fscanf(histdos_file, "%*s %lu", &rejectedMoves) != 1)
    std::cerr << "     ERROR! Cannot read rejectedMoves \n";

  if (fscanf(histdos_file, "%*s %d", &iterations) != 1)
    std::cerr << "     ERROR! Cannot read iterations \n";

  if (fscanf(histdos_file, "%*s %lu", &numBelowRange) != 1)
    std::cerr << "     ERROR! Cannot read numBelowRange \n";

  if (fscanf(histdos_file, "%*s %lu", &numAboveRange) != 1)
    std::cerr << "     ERROR! Cannot read numAboveRange \n";

  if (fscanf(histdos_file, "%*s %u", &numBinsFailingCriterion) != 1)
    std::cerr << "     ERROR! Cannot read numBinsFailingCriterion \n";

  if (fscanf(histdos_file, "%*s %d", &numHistogramNotImproved) != 1)
    std::cerr << "     ERROR! Cannot read numHistogramNotImproved \n";

  if (fscanf(histdos_file, "%*s %d", &numHistogramRefreshed) != 1)
    std::cerr << "     ERROR! Cannot read numHistogramRefreshed \n";

  // Open up arrays needed to store the mask, histogram and DOS
  hist.assign(numBins, 0);
  dos.assign(numBins, 0.0);
  visited.assign(numBins, 0);
  probDistribution.assign(numBins, 0);

  // Continue reading the histogram and DOS from file
  unsigned int dummy = 0;
  for (unsigned int i = 0; i < numBins; i++) {
    if (fscanf(histdos_file, "%u %d %lu %lf", &dummy, &visited[i], &hist[i], &dos[i]) != 4)
      std::cerr << "     ERROR! Cannot read histogram and DOS.\n";
    if (dummy != i)
      std::cerr << "     ERROR! Problem reading histogram and DOS. Check!\n";
  }

  //for (unsigned int i = 0; i < numBins; i++)
  //  printf("Check: %d %lu %20.8f\n", visited[i], hist[i], dos[i]);

  fclose(histdos_file);

/*
//YingWai's note:  (Sep 20, 16)
//The following works, but it gives a warning that dos[i] is an unsigned long int type... (!)

  std::ifstream HistDOSFile(fileName);
  std::string line, key;
  unsigned int i = 0;

  while (std::getline(HistDOSFile, line)) {
    if (!line.empty()) {
      std::istringstream lineStream(line);
      lineStream >> key;

      if (key.compare(0, 1, "#") != 0) {
        
        if (key == "dim") {
          lineStream >> dim;
          std::cerr <<  "YingWai's check for I/O. dim = " << dim << "\n";
          continue;
        }
        if (key == "flatnessCriterion") {
          lineStream >> flatnessCriterion;
          std::cerr <<  "YingWai's check for I/O. flatnessCriterion = " << flatnessCriterion << "\n";
          continue;
        }
        if (key == "modFactor") {
          lineStream >> modFactor;
          std::cerr <<  "YingWai's check for I/O. modFactor = " << modFactor << "\n";
          continue;
        }
        if (key == "modFactorFinal") {
          lineStream >> modFactorFinal;
          std::cerr <<  "YingWai's check for I/O. modFactorFinal = " << modFactorFinal << "\n";
          continue;
        }
        if (key == "modFactorReducer") {
          lineStream >> modFactorReducer;
          std::cerr <<  "YingWai's check for I/O. modFactorReducer = " << modFactorReducer << "\n";
          continue;
        }
        if (key == "histogramCheckInterval") {
          lineStream >> histogramCheckInterval;
          std::cerr <<  "YingWai's check for I/O. histogramCheckInterval = " << histogramCheckInterval << "\n";
          continue;
        }
        if (key == "histogramRefreshInterval") {
          lineStream >> histogramRefreshInterval;
          std::cerr <<  "YingWai's check for I/O. histogramRefreshInterval = " << histogramRefreshInterval << "\n";
          continue;
        }
        if (key == "Emin") {
          lineStream >> Emin;
          std::cerr <<  "YingWai's check for I/O. Emin = " << Emin << "\n";
          continue;
        }
        if (key == "Emax") {
          lineStream >> Emax;
          std::cerr <<  "YingWai's check for I/O. Emax = " << Emax << "\n";
          continue;
        }
        if (key == "binSize") {
          lineStream >> binSize;
          std::cerr <<  "YingWai's check for I/O. binSize = " << binSize << "\n";
          continue;
        }
        if (key == "numBins") {
          lineStream >> numBins;
          std::cerr <<  "YingWai's check for I/O. numBins = " << numBins << "\n";
          hist.assign(numBins, 0);
          dos.assign(numBins, 0.0);
          visited.assign(numBins, 0);
          continue;
        }
        if (key == "totalMCsteps") {
          lineStream >> totalMCsteps;
          std::cerr <<  "YingWai's check for I/O. totalMCsteps = " << totalMCsteps << "\n";
          continue;
        }
        if (key == "acceptedMoves") {
          lineStream >> acceptedMoves;
          std::cerr <<  "YingWai's check for I/O. acceptedMoves = " << acceptedMoves << "\n";
          continue;
        }
        if (key == "rejectedMoves") {
          lineStream >> rejectedMoves;
          std::cerr <<  "YingWai's check for I/O. rejectedMoves = " << rejectedMoves << "\n";
          continue;
        }
        if (key == "iterations") {
          lineStream >> iterations;
          std::cerr <<  "YingWai's check for I/O. iterations = " << iterations << "\n";
          continue;
        }
        if (key == "numBelowRange") {
          lineStream >> numBelowRange;
          std::cerr <<  "YingWai's check for I/O. numBelowRange = " << numBelowRange << "\n";
          continue;
        }
        if (key == "numAboveRange") {
          lineStream >> numAboveRange;
          std::cerr <<  "YingWai's check for I/O. numAboveRange = " << numAboveRange << "\n";
          continue;
        }
        if (key == "numBinsFailingCriterion") {
          lineStream >> numBinsFailingCriterion;
          std::cerr <<  "YingWai's check for I/O. numBinsFailingCriterion = " << numBinsFailingCriterion << "\n";
          continue;
        }
        // Continue reading these from file
        if (key == std::to_string(i)) {
          //std::cout << "key = " << stoul(key) << " = " << i << "\n";
          lineStream >> visited[i];
          lineStream >> hist[i];
          lineStream >> dos[i];
          //printf("YingWai's check for I/O : %d %d %20.5f\n", visited[i], hist[i], dos[i]);
          i++;
          //if (fscanf(histdos_file, "%u %d %lu %lf", &dummy, &visited[i], &hist[i], &dos[i]) != 4) 
          //  std::cerr << "Cannot read histogram and DOS.\n";
          continue;
        }

      }
    }
  }
*/

/*
  // YingWai's check, should be removed when things work fine
  std::cerr <<  "YingWai's check for I/O. dim = " << dim << "\n";
  std::cerr <<  "YingWai's check for I/O. flatnessCriterion = " << flatnessCriterion << "\n";
  std::cerr <<  "YingWai's check for I/O. modFactor = " << modFactor << "\n";
  std::cerr <<  "YingWai's check for I/O. modFactorFinal = " << modFactorFinal << "\n";
  std::cerr <<  "YingWai's check for I/O. modFactorReducer = " << modFactorReducer << "\n";
  std::cerr <<  "YingWai's check for I/O. histogramCheckInterval = " << histogramCheckInterval << "\n";
  std::cerr <<  "YingWai's check for I/O. Emin = " << Emin << "\n";
  std::cerr <<  "YingWai's check for I/O. Emax = " << Emax << "\n";
  std::cerr <<  "YingWai's check for I/O. binSize = " << binSize << "\n";
  std::cerr <<  "YingWai's check for I/O. numBins = " << numBins << "\n";
  std::cerr <<  "\n";

  std::cerr << "YingWai's check for I/O. totalMCsteps = " << totalMCsteps << "\n";
  std::cerr << "YingWai's check for I/O. acceptedMoves = " << acceptedMoves << "\n";
  std::cerr << "YingWai's check for I/O. rejectedMoves = " << rejectedMoves << "\n";
  std::cerr << "YingWai's check for I/O. iterations = " << iterations << "\n";
  std::cerr << "\n";

  // Write out histogram and DOS
  for (int i = 0; i < numBins; i++) {
    std::cerr << "YingWai's check for I/O. " 
              << i       << " " << visited[i] << " " 
              << hist[i] << " " << dos[i]     << "\n";
  }
*/

}


void Histogram::writeNormDOSFile(const char* fileName)
{

  FILE* dosFile;
  if (fileName != NULL) dosFile = fopen(fileName, "w");
  else dosFile = stdout;

  // To Do: see if there is a max function to use for C++ vector
  double maxDOS = 0.0;
  for (unsigned int i = 0; i < numBins; i++)
    if (dos[i] > maxDOS) maxDOS = dos[i];

  double norm = 0.0;
  for (unsigned int i = 0; i < numBins; i++)
    if (visited[i]) norm += exp(dos[i] - maxDOS);

  for (unsigned int i = 0; i < numBins; i++)
    fprintf(dosFile, "%18.10e  %18.10e\n", Emin + binSize * double(i),
                                           exp(dos[i] - maxDOS) / norm); 

  if (fileName != NULL) fclose(dosFile);

}


void Histogram::readMCInputFile(char const* fileName)
{

  if (GlobalComm.thisMPIrank == 0) 
    std::cout << "   Histogram class reading input file: " << fileName << "\n";

  std::ifstream inputFile(fileName);   // TODO: check if a file stream is initialized
  std::string line, key;

  if (inputFile.is_open()) {

    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
        
        std::istringstream lineStream(line);
        lineStream >> key;

        if (key.compare(0, 1, "#") != 0) {
          
          if (key == "dim") {
            lineStream >> dim;
            //std::cout << "WangLandau: dim = " << dim << "\n";
            continue;
          }
          if (key == "flatnessCriterion") {
            lineStream >> flatnessCriterion;
            //std::cout << "WangLandau: flatnessCriterion = " << flatnessCriterion << "\n";
            continue;
          }
          if (key == "modFactor") {
            lineStream >> modFactor;
            //std::cout << "WangLandau: modFactor = " << modFactor << "\n";
            continue;
          }
          if (key == "modFactorFinal") {
            lineStream >> modFactorFinal;
            //std::cout << "WangLandau: modFactorFinal = " << modFactorFinal << "\n";
            continue;
          }
          if (key == "modFactorReducer") {
            lineStream >> modFactorReducer;
            //std::cout << "WangLandau: modFactorReducer = " << modFactorReducer << "\n";
            continue;
          }
          if (key == "histogramCheckInterval") {
            lineStream >> histogramCheckInterval;
            //std::cout << "WangLandau: histogramCheckInterval = " << histogramCheckInterval << "\n";
            continue;
          }
          if (key == "histogramRefreshInterval") {
            lineStream >> histogramRefreshInterval;
            //std::cout << "WangLandau: histogramRefreshInterval = " << histogramRefreshInterval << "\n";
            continue;
          }
          if (key == "Emin") {
            lineStream >> Emin;
            //std::cout << "WangLandau: Emin = " << Emin << "\n";
            continue;
          }
          if (key == "Emax") {
            lineStream >> Emax;
            //std::cout << "WangLandau: Emax = " << Emax << "\n";
            continue;
          }
          if (key == "binSize") {
            lineStream >> binSize;
            //std::cout << "WangLandau: binSize = " << binSize << "\n";
            continue;
          }
          if (key == "numberOfWindows") {
            lineStream >> numberOfWindows;
            //std::cout << "REWL: numberOfWindows = " << numberOfWindows << "\n";
            continue;
          }
          if (key == "numberOfWalkersPerWindow") {
            lineStream >> numberOfWalkersPerWindow;
            //std::cout << "REWL: numberOfWalkersPerWindow = " << numberOfWalkersPerWindow << "\n";
            continue;
          }
          if (key == "overlap") {
            lineStream >> overlap;
            //std::cout << "REWL: overlap = " << overlap << "\n";
            continue;
          }
          if (key == "KullbackLeiblerDivergenceThreshold") {
            lineStream >> KullbackLeiblerDivergenceThreshold;
            //std::cout << "MUCA: KullbackLeiblerDivergenceThreshold = " << KullbackLeiblerDivergenceThreshold << "\n";
            continue;
          }
          if (key == "numberOfUpdatesPerIteration") {
            lineStream >> numberOfUpdatesPerIteration;
            //std::cout << "MUCA: numberOfUpdatesPerIteration = " << numberOfUpdatesPerIteration << "\n";
            continue;
          }
          if (key == "numberOfUpdatesMultiplier") {
            lineStream >> numberOfUpdatesMultiplier;
            //std::cout << "MUCA: numberOfUpdatesMultiplier = " << numberOfUpdatesMultiplier << "\n";
            continue;
          }
          if (key == "numberOfThermalizationSteps") {
            lineStream >> numberOfThermalizationSteps;
            //std::cout << "MUCA: numberOfThermalizationSteps = " << numberOfThermalizationSteps << "\n";
            continue;
          }

        }

      }

    }
    inputFile.close();

  }

}


