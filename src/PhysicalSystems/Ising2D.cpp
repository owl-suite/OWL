#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "Ising2D.hpp"
#include "Utilities/RandomNumberGenerator.hpp"

//Ising2D::Ising2D(SimulationInfo& sim_info, const char* filename, int initial)
Ising2D::Ising2D(const char* filename, int initial)
{

  printf("Simulation for 2D Ising model: %dx%d \n", simInfo.spinModelLatticeSize, simInfo.spinModelLatticeSize);

  int i, j;
  char c;

  Size = simInfo.spinModelLatticeSize;
  LatticeSize = Size * Size;
 
  spin = new SpinDirection[LatticeSize];

  //spin = new SpinDirection*[Size];
  //for (i = 0; i < Size; i++) 
  //  spin[i] = new SpinDirection[Size];

  if (filename != NULL) {
    FILE* f = fopen(filename, "r");
    if (f == NULL) {
      std::cout << "Coordinates file " << filename << " unreadable!" << std::endl;
      exit(1);
    }

    for(i = 0; i < Size; i++) {
      for (j = 0; j < Size; j++) {
        if (fscanf(f, "%c", &c) != 1) {
          std::cout << "Coordinates file " << filename << " unreadable!" << std::endl;
          exit(1);
        }
        switch (c) {
          case 'U' : { spin[i*Size+j] = 1; break; }
          default  : { spin[i*Size+j] = -1; }
          //case 'U' : { spin[i][j] = UP; break; }
          //default  : { spin[i][j] = DOWN; }
        }
      }
      fscanf(f, "%*c");
    }
    fclose(f);
  }
  else {
    for (i = 0; i < Size; i++) {
      for (j = 0; j < Size; j++) {

        switch (initial) {
        case 1  : {
          spin[i*Size+j] = -1;
          //spin[i][j] = DOWN;
	  break;
        }
        case 2  : {
          spin[i*Size+j] = 1;
          //spin[i][j] = UP;
	  break;
        }
        case 3  : {   // checkerboard
          if (((i + j) % 2) == 0) spin[i*Size+j] = -1;
          else spin[i*Size+j] = 1;
          //if (((i + j) % 2) == 0) spin[i][j] = DOWN;
          //else spin[i][j] = UP;
	  break;
        }
        default : {   // random
          if (getRandomNumber2() < 0.5) spin[i*Size+j] = -1;
          else spin[i*Size+j] = 1;
          //if (rng() < 0.5) spin[i][j] = DOWN;
          //else spin[i][j] = UP;
        }
        }

      }
    }
  }

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

  //for (int i = 0; i < Size; i++) 
  //  delete[] spin[i];
  delete[] spin;

  // Free MPI datatype
  pointerToConfiguration = NULL;
  MPI_Type_free(&MPI_ConfigurationType);

  deleteObservables();

  if (GlobalComm.thisMPIrank == 0)
    printf("Ising2D finished\n");

}



void Ising2D::writeConfiguration(int format, const char* filename)
{

  int x, y;

  FILE* f;
  if (filename != NULL) f = fopen(filename, "w");
  else f = stdout;

  switch (format) {

  default : {

    fprintf(f, "\n");
    fprintf(f, "2D Ising Model : %d x %d (%ld)\n", Size, Size, LatticeSize);
    fprintf(f, "Measures:");
    for (unsigned int i = 0; i < numObservables; i++)
      fprintf(f, " %10.5f", observables[i]);
      //fprintf(f, " %10d", observables[i]);
    fprintf(f, "\n");
    for (x = 0; x < Size; x++) {
      for (y = 0; y < Size; y++) 
        switch (spin[x*Size+y]) {
          //case UP : {fprintf(f, "U"); break;}
          case 1 : {fprintf(f, "U"); break;}
          default : {fprintf(f, "D");}
        }
      fprintf(f, "\n");
    }

  }

  }

  if (filename != NULL) fclose(f);

}


// YingWai: this does not seem to be needed anymore (Oct 10, 2017)
/*
void Ising2D::GetMeasuresBruteForce()
{
  //printf("!!! CALLING GetMeasuresBruteForce !!! \n");

  int x, y;
  int xLeft, yBelow;

  // Uncomment this when observables[] are used
  //resetObservables();
  int tempE = 0.0;
  int tempM = 0.0;

  for (x = 0; x < Size; x++) {
    if (x != 0) xLeft = x - 1; else xLeft = Size - 1;
    for (y = 0; y < Size; y++) {
      if (y != 0) yBelow = y - 1; else yBelow = Size - 1;
      tempE += spin[x][y] * (spin[xLeft][y] + spin[x][yBelow]);
      tempM += spin[x][y];
    }
  }
  tempE = -tempE;   // ferromagnetic interactions

  if (tempE != observables[0]) printf("Problem! tempE = %8d, observables[0] = %8.5f\n", tempE, observables[0]);
  if (tempM != observables[1]) printf("Problem! tempM = %8d, observables[1] = %8.5f\n", tempM, observables[1]);

}
*/


void Ising2D::getObservables()
{

  int x, y;
  int xLeft, yBelow;
  int xRight, yAbove;

  if (getObservablesFromScratch) {

    resetObservables();
  
    for (x = 0; x < Size; x++) {
      if (x != 0) xLeft = x - 1; else xLeft = Size - 1;
      for (y = 0; y < Size; y++) {
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
    int energyChange = sumNeighbor * CurType * 2;

    observables[0] += energyChange;
    observables[1] += spin[CurX*Size+CurY] - CurType;

    //printf("observables = %10.5f %10.5f \n", observables[0], observables[1]);
  }

  //GetMeasuresBruteForce();

}


void Ising2D::doMCMove()
{

  // Need this here since resetObservables() is not called if getObservablesFromScratch = false
  for (unsigned int i = 0; i < numObservables; i++)
    oldObservables[i] = observables[i];

  // randomly choose a site
  CurX = getIntRandomNumber() % Size;
  CurY = getIntRandomNumber() % Size;
  CurType = spin[CurX*Size + CurY];

  // flip the spin at that site
  if (CurType == -1)
    spin[CurX*Size+CurY] = 1;
  else
    spin[CurX*Size+CurY] = -1;
    
  //writeConfiguration(0);

}


/*
void Ising2D::undoMCMove()
{
  spin[CurX][CurY] = CurType;
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

  spin[CurX*Size+CurY] = CurType;
  for (unsigned int i=0; i < numObservables; i++)
    observables[i] = oldObservables[i];

}


void Ising2D::buildMPIConfigurationType()
{
 
  MPI_Type_contiguous(LatticeSize, MPI_INT, &MPI_ConfigurationType);
  MPI_Type_commit(&MPI_ConfigurationType);

}

