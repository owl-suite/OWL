#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "Ising2D.hpp"
#include "RandomNumberGenerator.hpp"

Ising2D::Ising2D(SimulationInfo& sim_info, const char* filename, int initial)
{

  printf("Simulation for 2D Ising model: %dx%d \n", sim_info.spinModelLatticeSize, sim_info.spinModelLatticeSize);

  int i, j;
  char c;

  Size = sim_info.spinModelLatticeSize;
  LatticeSize = Size * Size;

  spin = new SpinDirection*[Size];

  for (i = 0; i < Size; i++) 
    spin[i] = new SpinDirection[Size];

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
          case 'U' : { spin[i][j] = UP; break; }
          default  : { spin[i][j] = DOWN; }
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
          spin[i][j] = DOWN;
	  break;
        }
        case 2  : {
          spin[i][j] = UP;
	  break;
        }
        case 3  : {   // checkerboard
          if (((i + j) % 2) == 0) spin[i][j] = DOWN;
          else spin[i][j] = UP;
	  break;
        }
        default : {   // random
          if (rng() < 0.5) spin[i][j] = DOWN;
          else spin[i][j] = UP;
        }
        }

      }
    }
  }

  initializeObservables(2);      // observables[0] : energy
                                 // observables[1] : magnetization
  firstTimeGetMeasures = true;

}



Ising2D::~Ising2D()
{
  for (int i = 0; i < Size; i++) 
    delete[] spin[i];
  delete[] spin;

  deleteObservables();

  printf("Ising2D finished\n");
}


void Ising2D::readCommandLineOptions(SimulationInfo& sim_info)
{ };


void Ising2D::writeConfiguration(int format, const char* filename)
{

  int i;
  int x, y;

  FILE* f;
  if (filename != NULL) f = fopen(filename, "w");
  else f = stdout;

  switch (format) {

  default : {

    fprintf(f, "\n");
    fprintf(f, "2D Ising Model : %d x %d (%ld)\n", Size, Size, LatticeSize);
    fprintf(f, "Measures:");
    for (i = 0; i < numObservables; i++)
      fprintf(f, " %10d", observables[i]);
    fprintf(f, "\n");
    for (y = Size - 1; y >= 0; y--) {
      for (x = 0; x < Size; x++) 
        switch (spin[x][y]) {
          case UP : {fprintf(f, "U"); break;}
          default : {fprintf(f, "D");}
        }
      fprintf(f, "\n");
    }

  }

  }

  if (filename != NULL) fclose(f);

}



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

  if (tempE != observables[0]) printf("Problem! tempE = %8d, observables[0] = %8d\n", tempE, observables[0]);
  if (tempM != observables[1]) printf("Problem! tempM = %8d, observables[1] = %8d\n", tempM, observables[1]);

}



void Ising2D::getObservables()
{

  int x, y;
  int xLeft, yBelow;
  int xRight, yAbove;

  if (firstTimeGetMeasures) {

    resetObservables();
  
    for (x = 0; x < Size; x++) {
      if (x != 0) xLeft = x - 1; else xLeft = Size - 1;
      for (y = 0; y < Size; y++) {
        if (y != 0) yBelow = y - 1; else yBelow = Size - 1;
        observables[0] += spin[x][y] * (spin[xLeft][y] + spin[x][yBelow]);
        observables[1] += spin[x][y];
      }
    }
    observables[0] = -observables[0];     // ferromagnetic interaction
    firstTimeGetMeasures = false;
    printf("First time getObservables. \n");
  }
  else {
    if (CurX != 0) xLeft = CurX - 1; else xLeft = Size - 1;
    if (CurY != 0) yBelow = CurY - 1; else yBelow = Size - 1;
    if (CurX != (Size-1) ) xRight = CurX + 1; else xRight = 0;
    if (CurY != (Size-1) ) yAbove = CurY + 1; else yAbove = 0;

    int sumNeighbor = spin[xLeft][CurY] + spin[xRight][CurY] + spin[CurX][yBelow] + spin[CurX][yAbove];
    int energyChange = sumNeighbor * CurType * 2;

    observables[0] += energyChange;
    observables[1] += spin[CurX][CurY] - CurType;

    //printf("observables = %10.5f %10.5f \n", observables[0], observables[1]);
  }

  //GetMeasuresBruteForce();

}


void Ising2D::doMCMove()
{

  // Need this here since resetObservables() is not called if firstTimeGetMeasures = false
  for (int i = 0; i < numObservables; i++)
    oldObservables[i] = observables[i];

  CurX = rng() % Size;
  CurY = rng() % Size;

  CurType = spin[CurX][CurY];

  if (CurType == DOWN)
    spin[CurX][CurY] = UP;
  else
    spin[CurX][CurY] = DOWN;

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
  for (int i=0; i<numObservables; i++)
    oldObservables[i] = observables[i];
}


void Ising2D::rejectMCMove()
{
  spin[CurX][CurY] = CurType;
  for (int i=0; i<numObservables; i++)
    observables[i] = oldObservables[i];
}

