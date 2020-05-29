#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "CrystalStructure3D.hpp"
#include "Utilities/RandomNumberGenerator.hpp"


CrystalStructure3D::CrystalStructure3D(const char* inputFile, const char* spinConfigFile, int initial) : lattice(inputFile)
{

  printf("Simulation for customized 3D crystal structure: %dx%dx%d unit cells \n", lattice.unitCellDimensions[0], lattice.unitCellDimensions[1], lattice.unitCellDimensions[2]);

  assert (lattice.totalNumberOfAtoms > 0);
  spin = new SpinDirection[lattice.totalNumberOfAtoms];

  if (spinConfigFile != NULL)
    readSpinConfigFile(spinConfigFile);
  else
    initializeSpinConfiguration(initial);

  // ------OK up to here

  initializeObservables(4);      // observables[0] : total energy
                                 // observables[1] : magnetization in x-direction
                                 // observables[2] : magnetization in y-direction
                                 // observables[3] : magnetization in z-direction
  firstTimeGetMeasures = true;
  getObservables();

}


// OK
CrystalStructure3D::~CrystalStructure3D()
{

  delete spin;
  deleteObservables();

  printf("CrystalStructure3D finished\n");

}


//void CrystalStructure3D::readCommandLineOptions()
//{ };


// OK
void CrystalStructure3D::writeConfiguration(int format, const char* filename)
{

  FILE* f;
  if (filename != NULL) f = fopen(filename, "w");
  else f = stdout;

  switch (format) {

  default : {

    fprintf(f, "\n");
    fprintf(f, "Customized 3D crystal structure: %dx%dx%d unit cells \n", lattice.unitCellDimensions[0], lattice.unitCellDimensions[1], lattice.unitCellDimensions[2]);
    fprintf(f, "Total number of atoms: %u \n", lattice.totalNumberOfAtoms);
    fprintf(f, "Measures:");
    for (unsigned int i = 0; i < numObservables; i++)
      fprintf(f, " %10.5f", observables[i]);
    fprintf(f, "\n");

    for (unsigned int i = 0; i < lattice.totalNumberOfAtoms; i++)
          fprintf(f, "%8.5f %8.5f %8.5f\n", spin[i].x, spin[i].y, spin[i].z);

  }

  }

  if (filename != NULL) fclose(f);

}


// To implement
void CrystalStructure3D::GetMeasuresBruteForce() 
{

/*
  //printf("!!! CALLING GetMeasuresBruteForce !!! \n");

  int i, j, k;
  int xLeft, yBelow, zBackward;

  // Uncomment this when observables[] are used
  //resetObservables();

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

*/
}


// To implement
void CrystalStructure3D::getObservables() 
{


}

// OK
void CrystalStructure3D::doMCMove()
{

  // Need this here since resetObservables() is not called if firstTimeGetMeasures = false
  // for (int i = 0; i < numObservables; i++)
  //   oldObservables[i] = observables[i];

  currentPosition = getUnsignedIntRandomNumber() % lattice.totalNumberOfAtoms;
  currentSpin = spin[currentPosition];

  assignRandomSpinConfiguration(currentPosition);

}


/*
void CrystalStructure3D::undoMCMove()
{
  spin[CurX][CurY][CurZ] = currentSpin;
  restoreObservables();
}
*/

// OK
void CrystalStructure3D::acceptMCMove()
{
  // update "old" observables
  for (int i=0; i<numObservables; i++)
    oldObservables[i] = observables[i];
}

// OK
void CrystalStructure3D::rejectMCMove()
{
  spin[currentPosition] = currentSpin;
  for (int i=0; i<numObservables; i++)
    observables[i] = oldObservables[i];
}

/*
void CrystalStructure3D::buildMPIConfigurationType()
{
}
*/



void CrystalStructure3D::readSpinConfigFile(const char* spinConfigFile)
{

}



// OK
void CrystalStructure3D::initializeSpinConfiguration(int initial)
{

  for (unsigned int i = 0; i < lattice.totalNumberOfAtoms; i++) {

    switch (initial) {
      case 1 : {
        spin[i].x = 1.0;
        spin[i].y = 0.0;
        spin[i].z = 0.0;
        break;
      }
      case 2  : {
        spin[i].x = 0.0;
        spin[i].y = 1.0;
        spin[i].z = 0.0;
	    break;
      }
      case 3  : {
        spin[i].x = 0.0;
        spin[i].y = 0.0;
        spin[i].z = 1.0;
	      break;
      }
      case 4  : {
        spin[i].x = 0.0;
        spin[i].y = 0.0;
        spin[i].z = (i % 2 == 0) ? 1.0 : -1.0;
        break;
      }
      default  : {
        assignRandomSpinConfiguration(i);
      }
    }

  }

}

// OK
void CrystalStructure3D::assignRandomSpinConfiguration(unsigned int i)
{

  double r1, r2, rr;

  do {
    r1 = 2.0 * getRandomNumber();
    r2 = 2.0 * getRandomNumber();                
    rr = r1 * r1 + r2 * r2;
  } while (rr > 1.0);

  spin[i].x = 2.0 * r1 * sqrt(1.0 - rr);
  spin[i].y = 2.0 * r2 * sqrt(1.0 - rr);
  spin[i].z = 1.0 - 2.0 * rr;

}


// To implement
void CrystalStructure3D::readHamiltonianTerms(const char* inputFile)
{

  
}