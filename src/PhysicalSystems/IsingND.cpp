#include <cassert>
#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include "IsingND.hpp"
#include "Utilities/RandomNumberGenerator.hpp"


IsingND::IsingND(const char* spinConfigFile, int initial)
{

  printf("Simulation for %d-D Ising model of length %d \n", simInfo.spinModelDimension, simInfo.spinModelLatticeSize);

  Size      = simInfo.spinModelLatticeSize;
  dimension = simInfo.spinModelDimension;
  setSystemSize(Size, dimension);
  currentPosition.assign(dimension, 0);

  printf("YingWai's check: Size = %d, dimension = %d, systemSize = %d\n", Size, dimension, systemSize);

  spin = new IsingSpinDirection[systemSize];

  calculateOffsets();

  // Initialize configuration from file if applicable
  if (std::filesystem::exists(spinConfigFile))
    readSpinConfigFile(spinConfigFile);
  else if (simInfo.restartFlag && std::filesystem::exists("configurations/config_checkpoint.dat"))
    readSpinConfigFile("configurations/config_checkpoint.dat");
  else
    initializeSpinConfiguration(initial);

  observableName.push_back("Total energy, E");                            // observables[0] : total energy
  observableName.push_back("Total magnetization, M");                     // observables[1] : total magnetization
  observableName.push_back("Total absolute magnetization, |M|");          // observables[2] : total absolute magnetization
  initializeObservables(observableName.size());

  getObservablesFromScratch = true;
  getObservables();

  buildMPIConfigurationType();
  pointerToConfiguration = static_cast<void*>(&spin[0]);

}


IsingND::~IsingND()
{

  delete[] spin;

  // Free MPI datatype
  pointerToConfiguration = NULL;
  MPI_Type_free(&MPI_ConfigurationType);

  if (GlobalComm.thisMPIrank == 0)
    printf("\nIsingND finished\n");

}


void IsingND::writeConfiguration(int format, const char* filename)
{

  FILE* f;

  switch (format) {

    case 2 : {     // Write everything in one file, one line for each configuration

      if (filename != NULL) f = fopen(filename, "a");
      else f = stdout;

      // Write the configuration
      for (unsigned int i=0; i<systemSize; i++) {
        switch (spin[i]) {
          case 1  : {fprintf(f, "1 "); break;}
          default : {fprintf(f, "-1 ");}
        }
      }
      
      // Write the observables
      for (unsigned int i = 0; i < numObservables; i++)
        fprintf(f, " %10.5f", observables[i]);
      fprintf(f, "\n");

      break;
    }
    default : {

      if (filename != NULL) f = fopen(filename, "w");
      else f = stdout;

      fprintf(f, "# %u-D Ising Model of length %u \n\n", dimension, Size);
      fprintf(f, "TotalNumberOfSpins %u\n", systemSize);
      fprintf(f, "Observables ");
  
      for (unsigned int i = 0; i < numObservables; i++)
        fprintf(f, " %10.5f", observables[i]);
        //fprintf(f, "%50s %10.5f", observableName[i], observables[i]);
      fprintf(f, "\n");
   
      fprintf(f, "\nSpinConfiguration\n");
      for (unsigned int i=0; i<systemSize; i++) {
        switch (spin[i]) {
          case 1  : {fprintf(f, "1 "); break;}
          default : {fprintf(f, "-1 ");}
        }
      }
      fprintf(f, "\n");

    }
  }

  if (filename != NULL) fclose(f);

}


void IsingND::getObservables()
{

  // Neighbors' coordinates
  std::vector<unsigned int>  neighbor1Position (dimension, 0);
  std::vector<unsigned int>  neighbor2Position (dimension, 0);

  if (getObservablesFromScratch) {

    resetObservables();
    
    for (unsigned int i=0; i<systemSize; i++) {
      currentIndex = i;
      getCoordinatesFromIndex(currentIndex, currentPosition);
      neighbor1Position = currentPosition;

      IsingSpinDirection sumNeighbor = 0.0;
      for (unsigned int d=0; d<dimension; d++) {
        if (currentPosition[d] != 0) neighbor1Position[d] = currentPosition[d] - 1;
        else neighbor1Position[d] = Size - 1;
        sumNeighbor += spin[getIndexFromCoordinates(neighbor1Position)];
        neighbor1Position[d] = currentPosition[d];
      }
      observables[0] += ObservableType(spin[currentIndex] * sumNeighbor);
      observables[1] += ObservableType(spin[currentIndex]);
    }

    observables[0] = -observables[0];             // ferromagnetic interaction
    observables[2] = abs(observables[1]);
    getObservablesFromScratch = false;

  }
  else {

    neighbor1Position = currentPosition;
    neighbor2Position = currentPosition;
    IsingSpinDirection sumNeighbor = 0.0;

    for (unsigned int d=0; d<dimension; d++) {
      if (currentPosition[d] != 0) neighbor1Position[d] = currentPosition[d] - 1;
      else neighbor1Position[d] = Size - 1;
      if (currentPosition[d] != Size-1) neighbor2Position[d] = currentPosition[d] + 1;
      else neighbor2Position[d] = 0;
      sumNeighbor += spin[getIndexFromCoordinates(neighbor1Position)]
                   + spin[getIndexFromCoordinates(neighbor2Position)];
      neighbor1Position[d] = currentPosition[d];
      neighbor2Position[d] = currentPosition[d];
    }
    ObservableType energyChange = ObservableType(sumNeighbor * oldSpin * 2);

    observables[0] += energyChange;
    observables[1] += ObservableType(spin[currentIndex] - oldSpin);
    observables[2]  = abs(observables[1]);
    //printf("observables = %10.5f %10.5f %10.5f\n", observables[0], observables[1], observables[2]);
  }

}


void IsingND::doMCMove()
{

  // Need this here since resetObservables() is not called if getObservablesFromScratch = false
  for (unsigned int i = 0; i < numObservables; i++)
    oldObservables[i] = observables[i];

  // randomly choose a site
  currentIndex = unsigned(getIntRandomNumber()) % systemSize;
  getCoordinatesFromIndex(currentIndex, currentPosition);
  
  oldSpin = spin[currentIndex];

  // flip the spin at that site
  if (oldSpin == -1)
    spin[currentIndex] = 1;
  else
    spin[currentIndex] = -1;
    
  //writeConfiguration(0);

}


/*
void IsingND::undoMCMove()
{
  spin[CurX][CurY] = oldSpin;
  restoreObservables();
}
*/
//OK!
void IsingND::acceptMCMove()
{

  // update "old" observables
  for (unsigned int i=0; i < numObservables; i++)
    oldObservables[i] = observables[i];

}


void IsingND::rejectMCMove()
{

  spin[currentIndex] = oldSpin;
  for (unsigned int i=0; i < numObservables; i++)
    observables[i] = oldObservables[i];

}


void IsingND::buildMPIConfigurationType()
{
 
  MPI_Type_contiguous(int(systemSize), MPI_INT, &MPI_ConfigurationType);
  MPI_Type_commit(&MPI_ConfigurationType);

}


void IsingND::readSpinConfigFile(const std::filesystem::path& spinConfigFile)
{

  std::cout << "\n   IsingND class reading configuration file: " << spinConfigFile << "\n";

  std::ifstream inputFile(spinConfigFile);
  std::string line, key;
  unsigned int numberOfSpins {0};
  char c;

  if (inputFile.is_open()) {

    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
        std::istringstream lineStream(line);
        lineStream >> key;
        if (key.compare(0, 1, "#") != 0) {
          if (key == "TotalNumberOfSpins") {
            lineStream >> numberOfSpins;
            //std::cout << "   IsingND: numberOfSpins = " << numberOfSpins << "\n";
            continue;
          }
          else if (key == "Observables") {
            unsigned int counter = 0;
            while (lineStream && counter < numObservables) {
              lineStream >> observables[counter];
              //std::cout << "   IsingND: observables[" << counter << "] = " << observables[counter] << "\n";
              counter++;
            }
            continue;
          }
          else if (key == "SpinConfiguration") {
            //std::cout << "   IsingND: Spin Configuration read: \n";
            for (unsigned int i=0; i<systemSize; i++) {
              lineStream.clear();
              std::getline(inputFile, line);               
              if (!line.empty())  lineStream.str(line);
                lineStream >> c;
                switch (c) {
                  case 1  : { spin[i] = 1; break; }
                  default : { spin[i] = -1; }
                }
            }
            continue;
          }
        }

      }
    }

    inputFile.close();
  }

  // Sanity checks:
  assert(systemSize == numberOfSpins);

  printf("   Initial configuration read:\n");
  for (unsigned int i=0; i<Size; i++) {
    printf("   ");
    for (unsigned int j=0; j<Size; j++)
      printf("%3d ", spin[i*Size+j]);
    printf("\n");
  }

}


void IsingND::initializeSpinConfiguration(int initial)
{

  for (unsigned int i = 0; i < systemSize; i++) {

    switch (initial) {
      case 1  : {
        spin[i] = -1;
        break;
      }
      case 2  : {
        spin[i] = 1;
        break;
      }
      case 3  : {   // checkerboard
        if (((i + i/Size) % 2) == 0) spin[i] = -1;
        else spin[i] = 1;
        break;
      }
      default : {   // random
        if (getRandomNumber2() < 0.5) spin[i] = -1;
        else spin[i] = 1;
      }
    }

  }

}


indexType IsingND::getIndexFromCoordinates(unsigned int dim, ...)
{

  std::va_list coords;
  va_start(coords, dim);

  indexType index {0};
  for (unsigned int i = 0; i < dim; i++)
    index += va_arg(coords, int) * offsets[i];
  
  va_end(coords);
    
  return index;

}


indexType IsingND::getIndexFromCoordinates(std::vector<unsigned int> coords)
{
  
  assert(coords.size() == dimension);

  indexType index {0};
  for (unsigned int i = 0; i < dimension; i++) {
    index += coords[i] * offsets[i];
  }

  return index;

}


void IsingND::getCoordinatesFromIndex(indexType index, std::vector<unsigned int>& coords)
{

  indexType temp = index;
  for (int i = dimension-1; i >= 0; --i) {  
    coords[i] = temp / offsets[i];
    temp = temp % offsets[i];
  }

} 


void IsingND::calculateOffsets()
{

  offsets.assign(dimension, 1);
  
  for (unsigned int i = 1; i < dimension; i++)
    offsets[i] = offsets[i-1] * Size;

}