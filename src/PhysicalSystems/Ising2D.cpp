#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include "Ising2D.hpp"
#include "Utilities/RandomNumberGenerator.hpp"


Ising2D::Ising2D(const char* spinConfigFile, int initial)
{

  printf("Simulation for 2D Ising model: %dx%d \n", simInfo.spinModelLatticeSize, simInfo.spinModelLatticeSize);

  Size = simInfo.spinModelLatticeSize;
  setSystemSize(Size * Size);
 
  spin = new SpinDirection[systemSize];

  // Initialize configuration from file if applicable
  if (std::filesystem::exists(spinConfigFile))
    readSpinConfigFile(spinConfigFile);
  else if (simInfo.restartFlag && std::filesystem::exists("configurations/config_checkpoint.dat"))
    readSpinConfigFile("configurations/config_checkpoint.dat");
  else
    initializeSpinConfiguration(initial);

  initializeObservables(2);
  observableName.push_back("Total energy, E");                            // observables[0] : total energy
  observableName.push_back("Total magnetization, M");                     // observables[1] : total magnetization

  getObservablesFromScratch = true;
  getObservables();

  buildMPIConfigurationType();
  pointerToConfiguration = static_cast<void*>(&spin[0]);
 
}



Ising2D::~Ising2D()
{

  delete[] spin;

  // Free MPI datatype
  pointerToConfiguration = NULL;
  MPI_Type_free(&MPI_ConfigurationType);

  deleteObservables();

  if (GlobalComm.thisMPIrank == 0)
    printf("\nIsing2D finished\n");

}



void Ising2D::writeConfiguration(int format, const char* filename)
{

  FILE* f;
  if (filename != NULL) f = fopen(filename, "w");
  else f = stdout;

  switch (format) {

  default : {

    fprintf(f, "# 2D Ising Model : %u x %u\n\n", Size, Size);
    fprintf(f, "TotalNumberOfSpins %u)\n", systemSize);
    fprintf(f, "Observables ");

    for (unsigned int i = 0; i < numObservables; i++)
      fprintf(f, " %10.5f", observables[i]);
    fprintf(f, "\n");
 
    fprintf(f, "\nSpinConfiguration\n");
    for (unsigned int x = 0; x < Size; x++) {
      for (unsigned int y = 0; y < Size; y++) 
        switch (spin[x*Size+y]) {
          case 1  : {fprintf(f, "U"); break;}
          default : {fprintf(f, "D");}
        }
      fprintf(f, "\n");
    }

  }

  }

  if (filename != NULL) fclose(f);

}


void Ising2D::getObservables()
{

  unsigned int xLeft, yBelow;
  unsigned int xRight, yAbove;

  if (getObservablesFromScratch) {

    resetObservables();
  
    for (unsigned int x = 0; x < Size; x++) {
      if (x != 0) xLeft = x - 1; else xLeft = Size - 1;
      for (unsigned int y = 0; y < Size; y++) {
        if (y != 0) yBelow = y - 1; else yBelow = Size - 1;
        observables[0] += spin[x*Size+y] * (spin[xLeft*Size+y] + spin[x*Size+yBelow]);
        observables[1] += spin[x*Size+y];
      }
    }
    observables[0] = -observables[0];     // ferromagnetic interaction
    getObservablesFromScratch = false;
    //printf("Calculated observables from scratch. \n");
  }
  else {
    if (CurX != 0) xLeft = CurX - 1; else xLeft = Size - 1;
    if (CurY != 0) yBelow = CurY - 1; else yBelow = Size - 1;
    if (CurX != (Size-1) ) xRight = CurX + 1; else xRight = 0;
    if (CurY != (Size-1) ) yAbove = CurY + 1; else yAbove = 0;

    int sumNeighbor = spin[xLeft*Size+CurY] + spin[xRight*Size+CurY] + spin[CurX*Size+yBelow] + spin[CurX*Size+yAbove];
    int energyChange = sumNeighbor * oldSpin * 2;

    observables[0] += energyChange;
    observables[1] += spin[CurX*Size+CurY] - oldSpin;

    //printf("observables = %10.5f %10.5f \n", observables[0], observables[1]);
  }

}


void Ising2D::doMCMove()
{

  // Need this here since resetObservables() is not called if getObservablesFromScratch = false
  for (unsigned int i = 0; i < numObservables; i++)
    oldObservables[i] = observables[i];

  // randomly choose a site
  CurX = unsigned(getIntRandomNumber()) % Size;
  CurY = unsigned(getIntRandomNumber()) % Size;
  oldSpin = spin[CurX*Size + CurY];

  // flip the spin at that site
  if (oldSpin == -1)
    spin[CurX*Size+CurY] = 1;
  else
    spin[CurX*Size+CurY] = -1;
    
  //writeConfiguration(0);

}


/*
void Ising2D::undoMCMove()
{
  spin[CurX][CurY] = oldSpin;
  restoreObservables();
}
*/

void Ising2D::acceptMCMove()
{

  // update "old" observables
  for (unsigned int i=0; i < numObservables; i++)
    oldObservables[i] = observables[i];

}


void Ising2D::rejectMCMove()
{

  spin[CurX*Size+CurY] = oldSpin;
  for (unsigned int i=0; i < numObservables; i++)
    observables[i] = oldObservables[i];

}


void Ising2D::buildMPIConfigurationType()
{
 
  MPI_Type_contiguous(int(systemSize), MPI_INT, &MPI_ConfigurationType);
  MPI_Type_commit(&MPI_ConfigurationType);

}


void Ising2D::readSpinConfigFile(const std::filesystem::path& spinConfigFile)
{

  std::cout << "\n   Ising2D class reading configuration file: " << spinConfigFile << "\n";

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
            //std::cout << "   Ising2D: numberOfSpins = " << numberOfSpins << "\n";
            continue;
          }
          else if (key == "Observables") {
            unsigned int counter = 0;
            while (lineStream && counter < numObservables) {
              lineStream >> observables[counter];
              //std::cout << "   Ising2D: observables[" << counter << "] = " << observables[counter] << "\n";
              counter++;
            }
            continue;
          }
          else if (key == "SpinConfiguration") {
            //std::cout << "   Ising2D: Spin Configuration read: \n";
            for (unsigned int i=0; i<Size; i++) {
              lineStream.clear();
              std::getline(inputFile, line);               
              if (!line.empty())  lineStream.str(line);
              for (unsigned int j=0; j<Size; j++) {
                lineStream >> c;
                switch (c) {
                  case 'U' : { spin[i*Size+j] = 1; break; }
                  default  : { spin[i*Size+j] = -1; }
                }
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


void Ising2D::initializeSpinConfiguration(int initial)
{

  for (unsigned int i = 0; i < Size; i++) {
    for (unsigned int j = 0; j < Size; j++) {

      switch (initial) {
        case 1  : {
          spin[i*Size+j] = -1;
          break;
        }
        case 2  : {
          spin[i*Size+j] = 1;
          break;
        }
        case 3  : {   // checkerboard
          if (((i + j) % 2) == 0) spin[i*Size+j] = -1;
          else spin[i*Size+j] = 1;
          break;
        }
        default : {   // random
          if (getRandomNumber2() < 0.5) spin[i*Size+j] = -1;
          else spin[i*Size+j] = 1;
        }
      }

    }
  }

}