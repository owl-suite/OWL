#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include "Heisenberg3D.hpp"
#include "Utilities/RandomNumberGenerator.hpp"

Heisenberg3D::Heisenberg3D(const char* spinConfigFile, int initial)
{

  printf("Simulation for 3D Heisenberg model: %dx%dx%d \n", simInfo.spinModelLatticeSize, simInfo.spinModelLatticeSize, simInfo.spinModelLatticeSize);

  assert (simInfo.spinModelLatticeSize > 0);

  Size = simInfo.spinModelLatticeSize;
  setSystemSize(Size * Size * Size);

  spin = new SpinDirection**[Size];

  for (unsigned int i = 0; i < Size; i++) {
    spin[i] = new SpinDirection*[Size];
    for (unsigned int j = 0; j < Size; j++)
      spin[i][j] = new SpinDirection[Size];
  }

  if (std::filesystem::exists(spinConfigFile))
    readSpinConfigFile(spinConfigFile);
  else if (simInfo.restartFlag && std::filesystem::exists("configurations/config_checkpoint.dat"))
    readSpinConfigFile("configurations/config_checkpoint.dat");
  else
    initializeSpinConfiguration(initial);

  initializeObservables(6);
  observableName.push_back("Total energy, E");                            // observables[0] : total energy
  observableName.push_back("Magnetization in x-direction, M_x");          // observables[1] : magnetization in x-direction
  observableName.push_back("Magnetization in y-direction, M_y");          // observables[2] : magnetization in y-direction
  observableName.push_back("Magnetization in z-direction, M_z");          // observables[3] : magnetization in z-direction
  observableName.push_back("Total magnetization, M");                     // observables[4] : total magnetization
  observableName.push_back("4th order magnetization, M^4");               // observables[5] : total magnetization to the order 4

  firstTimeGetMeasures = true;
  getObservables();

}



Heisenberg3D::~Heisenberg3D()
{
  for (unsigned int i = 0; i < Size; i++) {
    for (unsigned int j = 0; j < Size; j++)
      delete[] spin[i][j];
    delete[] spin[i];
  }
  delete[] spin;

  deleteObservables();

  printf("Heisenberg3D finished\n");
}


//void Heisenberg3D::readCommandLineOptions()
//{ };


void Heisenberg3D::writeConfiguration(int format, const char* filename)
{

  FILE* f;
  if (filename != NULL) f = fopen(filename, "w");
  else f = stdout;

  switch (format) {

  default : {

    fprintf(f, "# 3D Heisenberg Model : %u x %u x %u \n\n", Size, Size, Size);
    fprintf(f, "TotalNumberOfSpins %u\n", systemSize);
    fprintf(f, "Observables ");

    for (unsigned int i = 0; i < numObservables; i++)
      fprintf(f, " %10.5f", observables[i]);
    fprintf(f, "\n");

    fprintf(f, "\nSpinConfiguration\n");
    for (unsigned int i = 0; i < Size; i++) {
      for (unsigned int j = 0; j < Size; j++) {
        for (unsigned int k = 0; k < Size; k++)
          fprintf(f, "%8.5f %8.5f %8.5f\n", spin[i][j][k].x, spin[i][j][k].y, spin[i][j][k].z);
      }
    }

  }

  }

  if (filename != NULL) fclose(f);

}


void Heisenberg3D::getObservables() 
{

  if (firstTimeGetMeasures) {
    //resetObservables();
    observables[0] = getExchangeInterations() + getExternalFieldEnergy();
    std::tie(observables[1], observables[2], observables[3], observables[4]) = getMagnetization();
    observables[5] = pow(observables[4], 4.0);

    firstTimeGetMeasures = false;
    //printf("First time getObservables. \n");
  }
  else {
    observables[0] += getDifferenceInExchangeInterations() + getDifferenceInExternalFieldEnergy();
    observables[1] += spin[CurX][CurY][CurZ].x - CurType.x;
    observables[2] += spin[CurX][CurY][CurZ].y - CurType.y;
    observables[3] += spin[CurX][CurY][CurZ].z - CurType.z;
    ObservableType temp = observables[1] * observables[1] + observables[2] * observables[2] + observables[3] * observables[3];
    observables[4] = sqrt(temp);
    observables[5] = temp * temp;
    //printf("observables = %10.5f %10.5f %10.5f %10.5f %10.5f\n", observables[0], observables[1], observables[2], observables[3], observables[4]);
  }

}


ObservableType Heisenberg3D::getExchangeInterations()
{
  
  unsigned int xLeft, yBelow, zBackward;
  ObservableType energy {0.0};

  for (unsigned int i = 0; i < Size; i++) {
    if (i != 0) xLeft = i - 1; else xLeft = Size - 1;
    for (unsigned int j = 0; j < Size; j++) {
      if (j != 0) yBelow = j - 1; else yBelow = Size - 1;
      for (unsigned int k = 0; k < Size; k++) {
        if (k != 0) zBackward = k - 1; else zBackward = Size - 1;

        energy += spin[i][j][k].x * (spin[xLeft][j][k].x + spin[i][yBelow][k].x + spin[i][j][zBackward].x) + 
                  spin[i][j][k].y * (spin[xLeft][j][k].y + spin[i][yBelow][k].y + spin[i][j][zBackward].y) +
                  spin[i][j][k].z * (spin[xLeft][j][k].z + spin[i][yBelow][k].z + spin[i][j][zBackward].z);

      }
    }  
  }
  
  return -energy; // ferromagnetic (FO) coupling

}


ObservableType Heisenberg3D::getExternalFieldEnergy()
{
  return 0.0;
}


std::tuple<ObservableType, ObservableType, ObservableType, ObservableType> Heisenberg3D::getMagnetization()
{

  ObservableType m1 {0.0};
  ObservableType m2 {0.0};
  ObservableType m3 {0.0};
  ObservableType m4 {0.0};

  for (unsigned int i = 0; i < Size; i++) {
    for (unsigned int j = 0; j < Size; j++) {
      for (unsigned int k = 0; k < Size; k++) {
        m1 += spin[i][j][k].x;
        m2 += spin[i][j][k].y;
        m3 += spin[i][j][k].z;
      }
    }
  }

  m4 = sqrt(m1 * m1 + m2 * m2 + m3 * m3);

  return {m1, m2, m3, m4};

}


ObservableType Heisenberg3D::getDifferenceInExchangeInterations()
{

  unsigned int xLeft, yBelow, zBackward;
  unsigned int xRight, yAbove, zForward;
  ObservableType energyChange {0.0};

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

  return -energyChange;           // ferromagnetic (FO) coupling

}


ObservableType Heisenberg3D::getDifferenceInExternalFieldEnergy()
{
  return 0.0;
}


void Heisenberg3D::doMCMove()
{

  double r1, r2, rr;

  // Need this here since resetObservables() is not called if firstTimeGetMeasures = false
  //for (int i = 0; i < numObservables; i++)
  //  oldObservables[i] = observables[i];

  CurX = unsigned(getIntRandomNumber()) % Size;
  CurY = unsigned(getIntRandomNumber()) % Size;
  CurZ = unsigned(getIntRandomNumber()) % Size;

  CurType = spin[CurX][CurY][CurZ];

  do {
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
  for (unsigned int i = 0; i < numObservables; i++)
    oldObservables[i] = observables[i];
}


void Heisenberg3D::rejectMCMove()
{
  spin[CurX][CurY][CurZ] = CurType;
  for (unsigned int i = 0; i < numObservables; i++)
    observables[i] = oldObservables[i];
}


/*
void Heisenberg3D::buildMPIConfigurationType()
{
}
*/


void Heisenberg3D::readSpinConfigFile(const std::filesystem::path& spinConfigFile)
{

  std::cout << "\n   Heisenberg3D class reading configuration file: " << spinConfigFile << "\n";

  std::ifstream inputFile(spinConfigFile);
  std::string line, key;
  unsigned int numberOfSpins {0};

  if (inputFile.is_open()) {

    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
        std::istringstream lineStream(line);
        lineStream >> key;

        if (key.compare(0, 1, "#") != 0) {

          if (key == "TotalNumberOfSpins") {
            lineStream >> numberOfSpins;
            //std::cout << "   Heisenberg3D: numberOfSpins = " << numberOfSpins << "\n";
            continue;
          }
          else if (key == "Observables") {
            unsigned int counter = 0;
            while (lineStream && counter < numObservables) {
              lineStream >> observables[counter];
              //std::cout << "   Heisenberg3D: observables[" << counter << "] = " << observables[counter] << "\n";
              counter++;
            }
            continue;
          }
          else if (key == "SpinConfiguration") {
            //std::cout << "   Heisenberg3D: Spin Configuration read: \n";
            for (unsigned int i=0; i<Size; i++) {
              for (unsigned int j=0; j<Size; j++) {
                for (unsigned int k=0; k<Size; k++) {
                  lineStream.clear();
                  std::getline(inputFile, line);               
                  if (!line.empty())  lineStream.str(line);
                  lineStream >> spin[i][j][k].x >> spin[i][j][k].y >> spin[i][j][k].z;
                  //printf("      %8.5f %8.5f %8.5f\n", spin[i][j][k].x, spin[i][j][k].y, spin[i][j][k].z);
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
  assert(numberOfSpins = systemSize);
  
  printf("   Initial configuration read:\n");
  for (unsigned int i=0; i<Size; i++) {
    for (unsigned int j=0; j<Size; j++) {
      for (unsigned int k=0; k<Size; k++)
        printf("      %8.5f %8.5f %8.5f\n", spin[i][j][k].x, spin[i][j][k].y, spin[i][j][k].z);
      printf("\n");
    }
    printf("\n");
  }

}


void Heisenberg3D::initializeSpinConfiguration(int initial)
{

  double r1, r2, rr;

  for (unsigned int i = 0; i < Size; i++) {
    for (unsigned int j = 0; j < Size; j++) {
      for (unsigned int k = 0; k < Size; k++) {

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
