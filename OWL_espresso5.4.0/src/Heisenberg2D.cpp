#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "Heisenberg2D.hpp"
#include "RandomNumberGenerator.hpp"

//Heisenberg2D::Heisenberg2D(SimulationInfo& sim_info, const char* filename, int initial)
Heisenberg2D::Heisenberg2D(const char* filename, int initial)
{

  printf("Simulation for 2D Heisenberg model: %dx%d \n", simInfo.spinModelLatticeSize, simInfo.spinModelLatticeSize);

  int i, j;
  char c;
  double r1, r2, rr;

  Size = simInfo.spinModelLatticeSize;
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
          case 'U' : {
            spin[i][j].x = 0.0;
            spin[i][j].y = 0.0;
            spin[i][j].z = 1.0;
            break;
          }
          default  : {
            spin[i][j].x = 0.0;
            spin[i][j].y = 0.0;
            spin[i][j].z = -1.0;
          }
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
          spin[i][j].x = 1.0;
          spin[i][j].y = 0.0;
          spin[i][j].z = 0.0;
	  break;
        }
        case 2  : {
          spin[i][j].x = 0.0;
          spin[i][j].y = 1.0;
          spin[i][j].z = 0.0;
	  break;
        }
        case 3  : {
          spin[i][j].x = 0.0;
          spin[i][j].y = 0.0;
          spin[i][j].z = 1.0;
	  break;
        }
        case 4  : {
          spin[i][j].x = 0.0;
          spin[i][j].y = 0.0;
	  if (((i + j) % 2) == 0) spin[i][j].z = 1.0;
	  else spin[i][j].z = -1.0;
	  break;
        }
        default : {
          do {
            //r1 = 1.0 - 2.0 * gsl_rng_uniform(rng);
            //r2 = 1.0 - 2.0 * gsl_rng_uniform(rng);
            r1 = 2.0 * getRandomNumber();
            r2 = 2.0 * getRandomNumber();
            rr = r1 * r1 + r2 * r2;
          } while (rr > 1.0);

          spin[i][j].x = 2.0 * r1 * sqrt(1.0 - rr);
          spin[i][j].y = 2.0 * r2 * sqrt(1.0 - rr);
          spin[i][j].z = 1.0 - 2.0 * rr;

        }
        }

      }
    }
  }

  initializeObservables(4);      // observables[0] : energy
                                 // observables[1] : magnetization in x-direction
                                 // observables[2] : magnetization in y-direction
                                 // observables[3] : magnetization in z-direction
  firstTimeGetMeasures = true;

}



Heisenberg2D::~Heisenberg2D()
{
  for (int i = 0; i < Size; i++) 
    delete[] spin[i];
  delete[] spin;

  deleteObservables();

  printf("Heisenberg2D finished\n");
}


void Heisenberg2D::readCommandLineOptions()
{ };


void Heisenberg2D::writeConfiguration(int format, const char* filename)
{

  int i;
  int x, y;

  FILE* f;
  if (filename != NULL) f = fopen(filename, "w");
  else f = stdout;

  switch (format) {

  default : {

    fprintf(f, "\n");
    fprintf(f, "2D Heisenberg Model : %d x %d (%ld)\n", Size, Size, LatticeSize);
    fprintf(f, "Measures:");
    for (i = 0; i < numObservables; i++)
      fprintf(f, " %10.5f", observables[i]);
    fprintf(f, "\n");
    for (y = Size - 1; y >= 0; y--) {
      for (x = 0; x < Size; x++) 
        fprintf(f, "%8.5f %8.5f %8.5f\n", spin[x][y].x, spin[x][y].y, spin[x][y].z);
      fprintf(f, "\n");
    }

  }

  }

  if (filename != NULL) fclose(f);

}



void Heisenberg2D::GetMeasuresBruteForce()
{
  //printf("!!! CALLING GetMeasuresBruteForce !!! \n");

  int x, y;
  int xLeft, yBelow;

  // Uncomment this when observables[] are used
  //resetObservables();
  double tempE = 0.0;
  double tempMx = 0.0;
  double tempMy = 0.0;
  double tempMz = 0.0;

  for (x = 0; x < Size; x++) {
    if (x != 0) xLeft = x - 1; else xLeft = Size - 1;
    for (y = 0; y < Size; y++) {
      if (y != 0) yBelow = y - 1; else yBelow = Size - 1;
      //observables[0] += spin[x][y].x * (spin[xLeft][y].x + spin[x][yBelow].x) + 
      //               spin[x][y].y * (spin[xLeft][y].y + spin[x][yBelow].y) +
      //               spin[x][y].z * (spin[xLeft][y].z + spin[x][yBelow].z);
      //observables[1] += spin[x][y].x;
      //observables[2] += spin[x][y].y;
      //observables[3] += spin[x][y].z;
      tempE  += spin[x][y].x * (spin[xLeft][y].x + spin[x][yBelow].x) + 
                spin[x][y].y * (spin[xLeft][y].y + spin[x][yBelow].y) +
                spin[x][y].z * (spin[xLeft][y].z + spin[x][yBelow].z);
      tempMx += spin[x][y].x;
      tempMy += spin[x][y].y;
      tempMz += spin[x][y].z;
    }
  }
  //observables[0] = -observables[0];   // ferromagnetic (FO) coupling
  tempE = -tempE;

  if ((abs(tempE) - abs(observables[0])) > 10e-8) printf("Problem! tempE - observables[0] = %15.10f\n", tempE-observables[0]);
  if ((abs(tempMx) - abs(observables[1])) > 10e-8) printf("Problem! tempMx - observables[1] = %15.10f\n", tempMx-observables[1]);
  if ((abs(tempMy) - abs(observables[2])) > 10e-8) printf("Problem! tempMy - observables[2] = %15.10f\n", tempMy-observables[2]);
  if ((abs(tempMz) - abs(observables[3])) > 10e-8) printf("Problem! tempMz - observables[3] = %15.10f\n", tempMz-observables[3]);

}



void Heisenberg2D::getObservables()
{

  int x, y;
  int xLeft, yBelow;
  int xRight, yAbove;
  double energyChange;

  if (firstTimeGetMeasures) {

    //resetObservables();
  
    for (x = 0; x < Size; x++) {
      if (x != 0) xLeft = x - 1; else xLeft = Size - 1;
      for (y = 0; y < Size; y++) {
        if (y != 0) yBelow = y - 1; else yBelow = Size - 1;
        observables[0] += spin[x][y].x * (spin[xLeft][y].x + spin[x][yBelow].x) + 
                          spin[x][y].y * (spin[xLeft][y].y + spin[x][yBelow].y) +
                          spin[x][y].z * (spin[xLeft][y].z + spin[x][yBelow].z);
        observables[1] += spin[x][y].x;
        observables[2] += spin[x][y].y;
        observables[3] += spin[x][y].z;
      }
    }
    observables[0] = -observables[0]; // ferromagnetic (FO) coupling
    firstTimeGetMeasures = false;
    printf("First time GetMeasures. \n");
  }
  else {
    if (CurX != 0) xLeft = CurX - 1; else xLeft = Size - 1;
    if (CurY != 0) yBelow = CurY - 1; else yBelow = Size - 1;
    if (CurX != (Size-1) ) xRight = CurX + 1; else xRight = 0;
    if (CurY != (Size-1) ) yAbove = CurY + 1; else yAbove = 0;

    energyChange = (spin[xLeft][CurY].x + spin[xRight][CurY].x + spin[CurX][yBelow].x + 
                    spin[CurX][yAbove].x) * (spin[CurX][CurY].x - CurType.x) +
                   (spin[xLeft][CurY].y + spin[xRight][CurY].y + spin[CurX][yBelow].y + 
                    spin[CurX][yAbove].y) * (spin[CurX][CurY].y - CurType.y) +
                   (spin[xLeft][CurY].z + spin[xRight][CurY].z + spin[CurX][yBelow].z + 
                    spin[CurX][yAbove].z) * (spin[CurX][CurY].z - CurType.z) ;
    energyChange = -energyChange;     // ferromagnetic (FO) coupling

    observables[0] += energyChange;
    observables[1] += spin[CurX][CurY].x - CurType.x;
    observables[2] += spin[CurX][CurY].y - CurType.y;
    observables[3] += spin[CurX][CurY].z - CurType.z;

    //printf("observables = %10.5f %10.5f %10.5f %10.5f\n", observables[0], observables[1], observables[2], observables[3]);
  }

  //GetMeasuresBruteForce();

}


void Heisenberg2D::doMCMove()
{

  double r1, r2, rr;

  // Need this here since resetObservables() is not called if firstTimeGetMeasures = false
  //for (int i = 0; i < numObservables; i++)
  //  oldObservables[i] = observables[i];

  //CurX = (int)(gsl_rng_uniform(rng) * Size);
  //CurY = (int)(gsl_rng_uniform(rng) * Size);
  CurX = rng() % Size;
  CurY = rng() % Size;

  CurType = spin[CurX][CurY];

  do {
    //r1 = 1.0 - 2.0 * gsl_rng_uniform(rng);
    //r2 = 1.0 - 2.0 * gsl_rng_uniform(rng);
    r1 = 2.0 * getRandomNumber();
    r2 = 2.0 * getRandomNumber();
    rr = r1 * r1 + r2 * r2;
  } while (rr > 1.0);

  spin[CurX][CurY].x = 2.0 * r1 * sqrt(1.0 - rr);
  spin[CurX][CurY].y = 2.0 * r2 * sqrt(1.0 - rr);
  spin[CurX][CurY].z = 1.0 - 2.0 * rr;

  //writeConfiguration(0);

}


/*
void Heisenberg2D::undoMCMove()
{
  spin[CurX][CurY] = CurType;
  restoreObservables();
}
*/

void Heisenberg2D::acceptMCMove()
{
  // update "old" observables
  for (int i=0; i<numObservables; i++)
    oldObservables[i] = observables[i];
}


void Heisenberg2D::rejectMCMove()
{
  spin[CurX][CurY] = CurType;
  for (int i=0; i<numObservables; i++)
    observables[i] = oldObservables[i];
}

