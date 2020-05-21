#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "Heisenberg3D.hpp"
#include "Utilities/RandomNumberGenerator.hpp"

Heisenberg3D::Heisenberg3D(const char* filename, int initial)
{

  printf("Simulation for 3D Heisenberg model: %dx%d \n", simInfo.spinModelLatticeSize, simInfo.spinModelLatticeSize);

  int i, j, k;
  double r1, r2, rr;

  Size = simInfo.spinModelLatticeSize;
  LatticeSize = Size * Size * Size;

  spin = new SpinDirection**[Size];

  for (i = 0; i < Size; i++) {
    spin[i] = new SpinDirection*[Size];
    for (j = 0; j < Size; j++)
      spin[i][j] = new SpinDirection[Size];
  }
  //SpinDirection spinTemp[LatticeSize];

  if (filename != NULL) {
    FILE* f = fopen(filename, "r");
    if (f == NULL) {
      std::cout << "Coordinates file " << filename << " unreadable!" << std::endl;
      exit(1);
    }

    for(i = 0; i < Size; i++) {
      for (j = 0; j < Size; j++) {
        for (k = 0; k < Size; k++) {
          //spinIndex = (long) i * Size * Size + j * Size + k;
          if (fscanf(f, "%lf %lf %lf", &spin[i][j][k].x, &spin[i][j][k].y, &spin[i][j][k].z) != 3) {
            std::cout << "Coordinates file " << filename << " unreadable!" << std::endl;
            exit(1);
          }
          //spin[i][j][k] = spinTemp[spinIndex];
          fscanf(f, "%*c");
        }
      }
    }
    fclose(f);
  }
  else {
    for (i = 0; i < Size; i++) {
      for (j = 0; j < Size; j++) {
        for (k = 0; k < Size; k++) {

          //spinIndex = (long) i * Size * Size + j * Size + k;
  
          switch (initial) {
            case 1 : {
              spin[i][j][k].x = 1.0;
              spin[i][j][k].y = 0.0;
              spin[i][j][k].z = 0.0;
              break;
            }
            case 2  : {
              spin[i][j][k].x = 0.0;
              spin[i][j][k].y = 1.0;
              spin[i][j][k].z = 0.0;
  	        break;
            }
            case 3  : {
              spin[i][j][k].x = 0.0;
              spin[i][j][k].y = 0.0;
              spin[i][j][k].z = 1.0;
  	          break;
            }
            case 4  : {
              spin[i][j][k].x = 0.0;
              spin[i][j][k].y = 0.0;
	            if (((i + j) % 2) == 0) spin[i][j][k].z = 1.0;
	            else spin[i][j][k].z = -1.0;
	            break;
            }
            default  : {
              do {
                //r1 = 1.0 - 2.0 * gsl_rng_uniform(rng);
                //r2 = 1.0 - 2.0 * gsl_rng_uniform(rng);
                r1 = 2.0 * getRandomNumber();
                r2 = 2.0 * getRandomNumber();                
                rr = r1 * r1 + r2 * r2;
              } while (rr > 1.0);

              spin[i][j][k].x = 2.0 * r1 * sqrt(1.0 - rr);
              spin[i][j][k].y = 2.0 * r2 * sqrt(1.0 - rr);
              spin[i][j][k].z = 1.0 - 2.0 * rr;
            }
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
  getObservables();

}



Heisenberg3D::~Heisenberg3D()
{
  for (int i = 0; i < Size; i++) {
    for (int j = 0; j < Size; j++)
      delete[] spin[i][j];
    delete[] spin[i];
  }
  delete[] spin;

  deleteObservables();

  printf("Heisenberg3D finished\n");
}


void Heisenberg3D::readCommandLineOptions()
{ };


void Heisenberg3D::writeConfiguration(int format, const char* filename)
{

  int i, j, k;

  FILE* f;
  if (filename != NULL) f = fopen(filename, "w");
  else f = stdout;

  switch (format) {

  default : {

    fprintf(f, "\n");
    fprintf(f, "3D Heisenberg Model : %d x %d x %d (%ld)\n", Size, Size, Size, LatticeSize);
    fprintf(f, "Measures:");
    for (i = 0; i < numObservables; i++)
      fprintf(f, " %10.5f", observables[i]);
    fprintf(f, "\n");
    for (i = 0; i < Size; i++) {
      for (j = 0; j < Size; j++) {
        for (k = 0; k < Size; k++)
          fprintf(f, "%8.5f %8.5f %8.5f\n", spin[i][j][k].x, spin[i][j][k].y, spin[i][j][k].z);
        fprintf(f, "\n");
      }
      fprintf(f, "\n");
    }

  }

  }

  if (filename != NULL) fclose(f);

}



void Heisenberg3D::GetMeasuresBruteForce() 
{
  //printf("!!! CALLING GetMeasuresBruteForce !!! \n");

  int i, j, k;
  int xLeft, yBelow, zBackward;

  // Uncomment this when observables[] are used
  //resetObservables();
  //double tempE = 0.0;
  //double tempMx = 0.0;
  //double tempMy = 0.0;
  //double tempMz = 0.0;
  ObservableType tempE = 0.0;
  ObservableType tempMx = 0.0;
  ObservableType tempMy = 0.0;
  ObservableType tempMz = 0.0;

  for (i = 0; i < Size; i++) {
    if (i != 0) xLeft = i - 1; else xLeft = Size - 1;
    for (j= 0; j < Size; j++) {
      if (j != 0) yBelow = j - 1; else yBelow = Size - 1;
      for (k = 0; k < Size; k++) {
        if (k != 0) zBackward = k - 1; else zBackward = Size - 1;
        //observables[0] += spin[x][y].x * (spin[xLeft][y].x + spin[x][yBelow].x) + 
        //               spin[x][y].y * (spin[xLeft][y].y + spin[x][yBelow].y) +
        //               spin[x][y].z * (spin[xLeft][y].z + spin[x][yBelow].z);
        //observables[1] += spin[x][y].x;
        //observables[2] += spin[x][y].y;
        //observables[3] += spin[x][y].z;
        tempE  += spin[i][j][k].x * (spin[xLeft][j][k].x + spin[i][yBelow][k].x + spin[i][j][zBackward].x) + 
                  spin[i][j][k].y * (spin[xLeft][j][k].y + spin[i][yBelow][k].y + spin[i][j][zBackward].y) +
                  spin[i][j][k].z * (spin[xLeft][j][k].z + spin[i][yBelow][k].z + spin[i][j][zBackward].z);
        tempMx += spin[i][j][k].x;
        tempMy += spin[i][j][k].y;
        tempMz += spin[i][j][k].z;
      }
    }
  }
  //observables[0] = -observables[0];   // ferromagnetic (FO) coupling
  tempE = -tempE;

  if ((std::abs(tempE) - std::abs(observables[0])) > 10e-8) printf("Problem! tempE - observables[0] = %15.10f\n", tempE-observables[0]);
  if ((std::abs(tempMx) - std::abs(observables[1])) > 10e-8) printf("Problem! tempMx - observables[1] = %15.10f\n", tempMx-observables[1]);
  if ((std::abs(tempMy) - std::abs(observables[2])) > 10e-8) printf("Problem! tempMy - observables[2] = %15.10f\n", tempMy-observables[2]);
  if ((std::abs(tempMz) - std::abs(observables[3])) > 10e-8) printf("Problem! tempMz - observables[3] = %15.10f\n", tempMz-observables[3]);

}



void Heisenberg3D::getObservables() 
{

  int i, j, k;
  int xLeft, yBelow, zBackward;
  int xRight, yAbove, zForward;
  //double energyChange;
  ObservableType energyChange;

  if (firstTimeGetMeasures) {

    //resetObservables();
  
    for (i = 0; i < Size; i++) {
      if (i != 0) xLeft = i - 1; else xLeft = Size - 1;
      for (j = 0; j < Size; j++) {
        if (j != 0) yBelow = j - 1; else yBelow = Size - 1;
        for (k = 0; k < Size; k++) {
          if (k != 0) zBackward = k - 1; else zBackward = Size - 1;
          observables[0] += spin[i][j][k].x * (spin[xLeft][j][k].x + spin[i][yBelow][k].x + spin[i][j][zBackward].x) + 
                            spin[i][j][k].y * (spin[xLeft][j][k].y + spin[i][yBelow][k].y + spin[i][j][zBackward].y) +
                            spin[i][j][k].z * (spin[xLeft][j][k].z + spin[i][yBelow][k].z + spin[i][j][zBackward].z);
          observables[1] += spin[i][j][k].x;
          observables[2] += spin[i][j][k].y;
          observables[3] += spin[i][j][k].z;
        }
      }  
    }
    observables[0] = -observables[0]; // ferromagnetic (FO) coupling
    firstTimeGetMeasures = false;
    printf("First time GetMeasures. \n");
  }
  else {
    if (CurX != 0) xLeft = CurX - 1; else xLeft = Size - 1;
    if (CurY != 0) yBelow = CurY - 1; else yBelow = Size - 1;
    if (CurZ != 0) zBackward = CurZ - 1; else zBackward = Size - 1;
    if (CurX != (Size-1) ) xRight = CurX + 1; else xRight = 0;
    if (CurY != (Size-1) ) yAbove = CurY + 1; else yAbove = 0;
    if (CurZ != (Size-1) ) zForward = CurZ + 1; else zForward = 0;

    energyChange = (spin[xLeft][CurY][CurZ].x + spin[xRight][CurY][CurZ].x + 
                    spin[CurX][yBelow][CurZ].x + spin[CurX][yAbove][CurZ].x + 
                    spin[CurX][CurY][zBackward].x + spin[CurX][CurY][zForward].x) * 
                   (spin[CurX][CurY][CurZ].x - CurType.x) +
                   (spin[xLeft][CurY][CurZ].y + spin[xRight][CurY][CurZ].y + 
                    spin[CurX][yBelow][CurZ].y + spin[CurX][yAbove][CurZ].y + 
                    spin[CurX][CurY][zBackward].y + spin[CurX][CurY][zForward].y) * 
                   (spin[CurX][CurY][CurZ].y - CurType.y) +
                   (spin[xLeft][CurY][CurZ].z + spin[xRight][CurY][CurZ].z + 
                    spin[CurX][yBelow][CurZ].z + spin[CurX][yAbove][CurZ].z + 
                    spin[CurX][CurY][zBackward].z + spin[CurX][CurY][zForward].z) * 
                   (spin[CurX][CurY][CurZ].z - CurType.z) ;
    energyChange = -energyChange;     // ferromagnetic (FO) coupling

    observables[0] += energyChange;
    observables[1] += spin[CurX][CurY][CurZ].x - CurType.x;
    observables[2] += spin[CurX][CurY][CurZ].y - CurType.y;
    observables[3] += spin[CurX][CurY][CurZ].z - CurType.z;

    //printf("observables = %10.5f %10.5f %10.5f %10.5f\n", observables[0], observables[1], observables[2], observables[3]);
  }

  //GetMeasuresBruteForce();

}


void Heisenberg3D::doMCMove()
{

  double r1, r2, rr;

  // Need this here since resetObservables() is not called if firstTimeGetMeasures = false
  //for (int i = 0; i < numObservables; i++)
  //  oldObservables[i] = observables[i];

  //CurX = (int)(gsl_rng_uniform(rng) * Size);
  //CurY = (int)(gsl_rng_uniform(rng) * Size);
  CurX = getIntRandomNumber() % Size;
  CurY = getIntRandomNumber() % Size;
  CurZ = getIntRandomNumber() % Size;

  CurType = spin[CurX][CurY][CurZ];

  do {
    //r1 = 1.0 - 2.0 * gsl_rng_uniform(rng);
    //r2 = 1.0 - 2.0 * gsl_rng_uniform(rng);
    r1 = 2.0 * getRandomNumber();
    r2 = 2.0 * getRandomNumber();
    rr = r1 * r1 + r2 * r2;
  } while (rr > 1.0);

  spin[CurX][CurY][CurZ].x = 2.0 * r1 * sqrt(1.0 - rr);
  spin[CurX][CurY][CurZ].y = 2.0 * r2 * sqrt(1.0 - rr);
  spin[CurX][CurY][CurZ].z = 1.0 - 2.0 * rr;

  //writeConfiguration(0);

}


/*
void Heisenberg3D::undoMCMove()
{
  spin[CurX][CurY][CurZ] = CurType;
  restoreObservables();
}
*/

void Heisenberg3D::acceptMCMove()
{
  // update "old" observables
  for (int i=0; i<numObservables; i++)
    oldObservables[i] = observables[i];
}


void Heisenberg3D::rejectMCMove()
{
  spin[CurX][CurY][CurZ] = CurType;
  for (int i=0; i<numObservables; i++)
    observables[i] = oldObservables[i];
}


void Heisenberg3D::buildMPIConfigurationType()
{
}

