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
//Histogram::Histogram(int restart, const char* inputFile)
Histogram::Histogram(int restart, const char* inputFile, const char* checkPointFile)
{

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
    std::cout << "Error: No input file for reading histogram's info. Quiting... \n";
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
    numBins  = ceil((Emax - Emin) / binSize) + 1;

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


int Histogram::getNumberOfBins()
{
  return numBins;
}


double Histogram::getDOS(ObservableType energy)
{
  idx = getIndex(energy);
  if (idx >= 0)
    return dos[idx];
  else {
    std::cerr << "Problem in File " << __FILE__ << " Line " << __LINE__  << std::endl;
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


void Histogram::setNumberOfBins(long int n)
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
    // If it is the first time a bin is visited:
    //   1. see if it can reference the DOS from neighboring bins
    //   2. reset Histogram and start over
    if ( visited[idx] == 0 ) {
      int refIdx = idx;
      if ( static_cast<unsigned>(idx) == 0 )
        refIdx = 1;
      else if ( static_cast<unsigned>(idx) == (numBins-1) )
        refIdx = idx - 1;
      else {
        if ( (visited[idx-1] > 0) && (visited[idx+1] > 0) )
          refIdx = ( dos[idx-1] < dos[idx+1] ? (idx-1) : (idx+1) );
        else if ( visited[idx-1] == 0 )
          refIdx = idx + 1;
        else if ( visited[idx+1] == 0 )
          refIdx = idx - 1;
      }

      dos[idx] = dos[refIdx];
      visited[idx] = 1;
      refreshHistogram();
    }
    else {
      dos[idx] += modFactor;
      hist[idx]++;
      //visited[idx] = 1;
    }
  }
  else { 
    std::cerr << "Error: idx < 0 in updateHistogramDOS!!\n";
    std::cerr << "Aborting...\n";
    exit(10);
  }
  //std::cerr << "energy = " << energy << std::endl;
  //std::cerr << "idx = " << idx << std::endl;
  //std::cerr << "visited[idx] = " << visited[idx] << std::endl;
  //std::cerr << "hist[idx] = " << hist[idx] << std::endl;
}

void Histogram::updateHistogram(ObservableType energy)
{
  idx = getIndex(energy);
  hist[idx]++;
  visited[idx] = 1;
  //std::cerr << "idx = " << idx << std::endl;
  //std::cerr << "visited[idx] = " << visited[idx] << std::endl;
}

void Histogram::updateDOS(ObservableType energy)
{
  idx = getIndex(energy);
  dos[idx] += modFactor;
  visited[idx] = 1;
  //std::cerr << "idx = " << idx << std::endl;
  //std::cerr << "visited[idx] = " << visited[idx] << std::endl;
}

void Histogram::updateDOSwithHistogram()
{
  for (unsigned int i=0; i<numBins; i++)
    if ((visited[i] == 1) && (hist[i] != 0))
      dos[i] += log(hist[i]);
}

bool Histogram::checkEnergyInRange(ObservableType energy)
{
  bool isWithinRange {false};
  if (energy < Emin) {
    //std::cerr << "Energy below range. Energy = " << energy << std::endl;
    numBelowRange++; 
    isWithinRange = false;
  }
  else if (energy > Emax) {
    //std::cerr << "Energy above range. Energy = " << energy << std::endl;
    numAboveRange++;
    isWithinRange = false;
  }
  else {
    //std::cerr << "Energy within range. Energy = " << energy << std::endl;
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
      if (hist[i] < flatnessReference) {
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
  int numVisitedBins = std::count(visited.begin(), visited.end(), 1);
  //int numVisitedBins = 0;
  //for (unsigned int i=0; i<numBins; i++) {
  //  if (visited[i] == 1)
  //    numVisitedBins++;
  //}
  std::cout << "Number of visited bins = " << numVisitedBins << std::endl;

  double flatnessReference = 1.0 / static_cast<double>( numVisitedBins );
  //double flatnessReference = 1.0 / static_cast<double>( std::max(numVisitedBins, 10) );

  KullbackLeiblerDivergence = 0.0;
  for (unsigned int i=0; i<numBins; i++) {
    if (visited[i] == 1) {
      probDistribution[i] = static_cast<double>(hist[i]) / static_cast<double>(numberOfUpdatesPerIteration);
      KullbackLeiblerDivergence += probDistribution[i] * log(probDistribution[i]/flatnessReference);
    }
  }

  std::cout << "KullbackLeiblerDivergence = " << KullbackLeiblerDivergence << std::endl;

  if (KullbackLeiblerDivergence <= KullbackLeiblerDivergenceThreshold)
    return true;
  else  
    return false;

}


// TO DO: construct filename from iteration and walkerID  (Sep 25, 2017)
void Histogram::writeHistogramDOSFile(const char* fileName, int iteration, int walkerID)
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
    //std::cerr << "visited[" << i << "] = " << visited[i] << std::endl;
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
  return floor(static_cast<double>(energy - Emin) / static_cast<double>(binSize));
}


void Histogram::readHistogramDOSFile(const char* fileName)
{

  if (GlobalComm.thisMPIrank == 0)
    std::cout << "Reading histogram checkpoint file : " << fileName << std::endl;

  FILE *histdos_file = fopen(fileName, "r");
  if (histdos_file == NULL) {
    std::cerr << "Error: cannot open histogram checkpoint file "  << fileName << std::endl;
    exit(7);    // perhaps can start from wl.input instead of quitting?
  }

  // TO DO: need to error-proof if they are not in order / missing...
  if (fscanf(histdos_file, "%*s %d", &dim) != 1)
    std::cerr << "Cannot read dim \n";    

  if (fscanf(histdos_file, "%*s %lf", &flatnessCriterion) != 1)
    std::cerr << "Cannot read flatnessCriterion \n";

  if (fscanf(histdos_file, "%*s %lf", &modFactor) != 1)
    std::cerr << "Cannot read modFactor \n";

  double modFactorFinal_tmp;
  if (fscanf(histdos_file, "%*s %lf", &modFactorFinal_tmp) != 1)
    std::cerr << "Cannot read modFactorFinal \n";
  if (modFactorFinal_tmp < modFactorFinal)
    std::cout << "modFactorFinal read from the checkpoint file is smaller than the input file. Simulation will continue with the new modFactorFinal. \n";

  if (fscanf(histdos_file, "%*s %lf", &modFactorReducer) != 1)
    std::cerr << "Cannot read modFactorReducer \n";

  if (fscanf(histdos_file, "%*s %u", &histogramCheckInterval) != 1)
    std::cerr << "Cannot read histogramCheckInterval \n";

  if (fscanf(histdos_file, "%*s %u", &histogramRefreshInterval) != 1)
    std::cerr << "Cannot read histogramRefreshInterval \n";

  if (fscanf(histdos_file, "%*s %lf", &Emin) != 1)
    std::cerr << "Cannot read Emin \n";

  if (fscanf(histdos_file, "%*s %lf", &Emax) != 1)
    std::cerr << "Cannot read Emax \n";

  if (fscanf(histdos_file, "%*s %lf", &binSize) != 1)
    std::cerr << "Cannot read binSize \n";

  if (fscanf(histdos_file, "%*s %u", &numBins) != 1)
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

  if (fscanf(histdos_file, "%*s %u", &numBinsFailingCriterion) != 1)
    std::cerr << "Cannot read numBinsFailingCriterion \n";

  if (fscanf(histdos_file, "%*s %d", &numHistogramNotImproved) != 1)
    std::cerr << " cannot read numHistogramNotImproved \n";

  if (fscanf(histdos_file, "%*s %d", &numHistogramRefreshed) != 1)
    std::cerr << " cannot read numHistogramRefreshed \n";

  // Open up arrays needed to store the mask, histogram and DOS
  hist.assign(numBins, 0);
  dos.assign(numBins, 0.0);
  visited.assign(numBins, 0);
  probDistribution.assign(numBins, 0);

  // Continue reading the histogram and DOS from file
  unsigned int dummy = 0;
  for (unsigned int i = 0; i < numBins; i++) {
    if (fscanf(histdos_file, "%u %d %lu %lf", &dummy, &visited[i], &hist[i], &dos[i]) != 4)
      std::cerr << "Cannot read histogram and DOS.\n";
    if (dummy != i)
      std::cerr << "Problem reading histogram and DOS. Check!\n";
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
          std::cerr <<  "YingWai's check for I/O. dim = " << dim << std::endl;
          continue;
        }
        if (key == "flatnessCriterion") {
          lineStream >> flatnessCriterion;
          std::cerr <<  "YingWai's check for I/O. flatnessCriterion = " << flatnessCriterion << std::endl;
          continue;
        }
        if (key == "modFactor") {
          lineStream >> modFactor;
          std::cerr <<  "YingWai's check for I/O. modFactor = " << modFactor << std::endl;
          continue;
        }
        if (key == "modFactorFinal") {
          lineStream >> modFactorFinal;
          std::cerr <<  "YingWai's check for I/O. modFactorFinal = " << modFactorFinal << std::endl;
          continue;
        }
        if (key == "modFactorReducer") {
          lineStream >> modFactorReducer;
          std::cerr <<  "YingWai's check for I/O. modFactorReducer = " << modFactorReducer << std::endl;
          continue;
        }
        if (key == "histogramCheckInterval") {
          lineStream >> histogramCheckInterval;
          std::cerr <<  "YingWai's check for I/O. histogramCheckInterval = " << histogramCheckInterval << std::endl;
          continue;
        }
        if (key == "histogramRefreshInterval") {
          lineStream >> histogramRefreshInterval;
          std::cerr <<  "YingWai's check for I/O. histogramRefreshInterval = " << histogramRefreshInterval << std::endl;
          continue;
        }
        if (key == "Emin") {
          lineStream >> Emin;
          std::cerr <<  "YingWai's check for I/O. Emin = " << Emin << std::endl;
          continue;
        }
        if (key == "Emax") {
          lineStream >> Emax;
          std::cerr <<  "YingWai's check for I/O. Emax = " << Emax << std::endl;
          continue;
        }
        if (key == "binSize") {
          lineStream >> binSize;
          std::cerr <<  "YingWai's check for I/O. binSize = " << binSize << std::endl;
          continue;
        }
        if (key == "numBins") {
          lineStream >> numBins;
          std::cerr <<  "YingWai's check for I/O. numBins = " << numBins << std::endl;
          hist.assign(numBins, 0);
          dos.assign(numBins, 0.0);
          visited.assign(numBins, 0);
          continue;
        }
        if (key == "totalMCsteps") {
          lineStream >> totalMCsteps;
          std::cerr <<  "YingWai's check for I/O. totalMCsteps = " << totalMCsteps << std::endl;
          continue;
        }
        if (key == "acceptedMoves") {
          lineStream >> acceptedMoves;
          std::cerr <<  "YingWai's check for I/O. acceptedMoves = " << acceptedMoves << std::endl;
          continue;
        }
        if (key == "rejectedMoves") {
          lineStream >> rejectedMoves;
          std::cerr <<  "YingWai's check for I/O. rejectedMoves = " << rejectedMoves << std::endl;
          continue;
        }
        if (key == "iterations") {
          lineStream >> iterations;
          std::cerr <<  "YingWai's check for I/O. iterations = " << iterations << std::endl;
          continue;
        }
        if (key == "numBelowRange") {
          lineStream >> numBelowRange;
          std::cerr <<  "YingWai's check for I/O. numBelowRange = " << numBelowRange << std::endl;
          continue;
        }
        if (key == "numAboveRange") {
          lineStream >> numAboveRange;
          std::cerr <<  "YingWai's check for I/O. numAboveRange = " << numAboveRange << std::endl;
          continue;
        }
        if (key == "numBinsFailingCriterion") {
          lineStream >> numBinsFailingCriterion;
          std::cerr <<  "YingWai's check for I/O. numBinsFailingCriterion = " << numBinsFailingCriterion << std::endl;
          continue;
        }
        // Continue reading these from file
        if (key == std::to_string(i)) {
          //std::cout << "key = " << stoul(key) << " = " << i << std::endl;
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


// TO DO: construct filename from walkerID  (Sep 25, 2017)
void Histogram::writeNormDOSFile(const char* fileName, int walkerID)
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
    std::cout << "Histogram class reading input file: " << fileName << std::endl;

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
            //std::cout << "WangLandau: dim = " << dim << std::endl;
            continue;
          }
          if (key == "flatnessCriterion") {
            lineStream >> flatnessCriterion;
            //std::cout << "WangLandau: flatnessCriterion = " << flatnessCriterion << std::endl;
            continue;
          }
          if (key == "modFactor") {
            lineStream >> modFactor;
            //std::cout << "WangLandau: modFactor = " << modFactor << std::endl;
            continue;
          }
          if (key == "modFactorFinal") {
            lineStream >> modFactorFinal;
            //std::cout << "WangLandau: modFactorFinal = " << modFactorFinal << std::endl;
            continue;
          }
          if (key == "modFactorReducer") {
            lineStream >> modFactorReducer;
            //std::cout << "WangLandau: modFactorReducer = " << modFactorReducer << std::endl;
            continue;
          }
          if (key == "histogramCheckInterval") {
            lineStream >> histogramCheckInterval;
            //std::cout << "WangLandau: histogramCheckInterval = " << histogramCheckInterval << std::endl;
            continue;
          }
          if (key == "histogramRefreshInterval") {
            lineStream >> histogramRefreshInterval;
            //std::cout << "WangLandau: histogramRefreshInterval = " << histogramRefreshInterval << std::endl;
            continue;
          }
          if (key == "Emin") {
            lineStream >> Emin;
            //std::cout << "WangLandau: Emin = " << Emin << std::endl;
            continue;
          }
          if (key == "Emax") {
            lineStream >> Emax;
            //std::cout << "WangLandau: Emax = " << Emax << std::endl;
            continue;
          }
          if (key == "binSize") {
            lineStream >> binSize;
            //std::cout << "WangLandau: binSize = " << binSize << std::endl;
            continue;
          }
          if (key == "numberOfWindows") {
            lineStream >> numberOfWindows;
            //std::cout << "REWL: numberOfWindows = " << numberOfWindows << std::endl;
            continue;
          }
          if (key == "numberOfWalkersPerWindow") {
            lineStream >> numberOfWalkersPerWindow;
            //std::cout << "REWL: numberOfWalkersPerWindow = " << numberOfWalkersPerWindow << std::endl;
            continue;
          }
          if (key == "overlap") {
            lineStream >> overlap;
            //std::cout << "REWL: overlap = " << overlap << std::endl;
            continue;
          }
          if (key == "KullbackLeiblerDivergenceThreshold") {
            lineStream >> KullbackLeiblerDivergenceThreshold;
            //std::cout << "MUCA: KullbackLeiblerDivergenceThreshold = " << KullbackLeiblerDivergenceThreshold << std::endl;
            continue;
          }
          if (key == "numberOfUpdatesPerIteration") {
            lineStream >> numberOfUpdatesPerIteration;
            //std::cout << "MUCA: numberOfUpdatesPerIteration = " << numberOfUpdatesPerIteration << std::endl;
            continue;
          }
          if (key == "numberOfUpdatesMultiplier") {
            lineStream >> numberOfUpdatesMultiplier;
            //std::cout << "MUCA: numberOfUpdatesMultiplier = " << numberOfUpdatesMultiplier << std::endl;
            continue;
          }
          if (key == "numberOfThermalizationSteps") {
            lineStream >> numberOfThermalizationSteps;
            //std::cout << "MUCA: numberOfThermalizationSteps = " << numberOfThermalizationSteps << std::endl;
            continue;
          }

        }

      }

    }
    inputFile.close();

  }

}


