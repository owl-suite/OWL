#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include "Heisenberg2D.hpp"
#include "Utilities/RandomNumberGenerator.hpp"

Heisenberg2D::Heisenberg2D(const char* spinConfigFile, int initial)
{

  printf("Simulation for 2D Heisenberg model: %dx%d \n", simInfo.spinModelLatticeSize, simInfo.spinModelLatticeSize);

  Size = simInfo.spinModelLatticeSize;
  setSystemSize(Size * Size);

  spin = new SpinDirection*[Size];

  for (unsigned int i = 0; i < Size; i++) 
    spin[i] = new SpinDirection[Size];

  if (std::filesystem::exists(spinConfigFile))
    readSpinConfigFile(spinConfigFile);
  else if (simInfo.restartFlag && std::filesystem::exists("configurations/config_checkpoint.dat"))
    readSpinConfigFile("configurations/config_checkpoint.dat");
  else
    initializeSpinConfiguration(initial);

  initializeObservables(5);
  observableName.push_back("Total energy, E");                            // observables[0] : total energy
  observableName.push_back("Magnetization in x-direction, M_x");          // observables[1] : magnetization in x-direction
  observableName.push_back("Magnetization in y-direction, M_y");          // observables[2] : magnetization in y-direction
  observableName.push_back("Magnetization in z-direction, M_z");          // observables[3] : magnetization in z-direction
  observableName.push_back("Total magnetization, M");                     // observables[4] : total magnetization

  firstTimeGetMeasures = true;
  getObservables();

}



Heisenberg2D::~Heisenberg2D()
{
  for (unsigned int i = 0; i < Size; i++) 
    delete[] spin[i];
  delete[] spin;

  deleteObservables();

  printf("Heisenberg2D finished\n");
}


//void Heisenberg2D::readCommandLineOptions()
//{ };


void Heisenberg2D::writeConfiguration(int format, const char* filename)
{

  FILE* f;
  if (filename != NULL) f = fopen(filename, "w");
  else f = stdout;

  switch (format) {

  default : {

    fprintf(f, "# 2D Heisenberg Model : %u x %u \n\n", Size, Size);
    fprintf(f, "TotalNumberOfSpins %u\n", systemSize);
    fprintf(f, "Observables ");

    for (unsigned int i = 0; i < numObservables; i++)
      fprintf(f, " %10.5f", observables[i]);
    fprintf(f, "\n");

    fprintf(f, "\nSpinConfiguration\n");
    for (unsigned int i = 0; i < Size; i++) {
      for (unsigned int j = 0; j < Size; j++)
        fprintf(f, "%8.5f %8.5f %8.5f\n", spin[i][j].x, spin[i][j].y, spin[i][j].z);
    }

  }

  }

  if (filename != NULL) fclose(f);

}


void Heisenberg2D::getObservables()
{

  if (firstTimeGetMeasures) {
    //resetObservables();
    observables[0] = getExchangeInteractions() + getExternalFieldEnergy();
    std::tie(observables[1], observables[2], observables[3], observables[4]) = getMagnetization();

    firstTimeGetMeasures = false;
    //printf("First time getObservables. \n");
  }
  else {
    observables[0] += getDifferenceInExchangeInteractions() + getDifferenceInExternalFieldEnergy();
    observables[1] += spin[CurX][CurY].x - CurType.x;
    observables[2] += spin[CurX][CurY].y - CurType.y;
    observables[3] += spin[CurX][CurY].z - CurType.z;
    observables[4] = sqrt(observables[1] * observables[1] + observables[2] * observables[2] + observables[3] * observables[3]);

    //printf("observables = %10.5f %10.5f %10.5f %10.5f %10.5f\n", observables[0], observables[1], observables[2], observables[3], observables[4]);
  }

}


ObservableType Heisenberg2D::getExchangeInteractions()
{

  unsigned int xLeft, yBelow;
  ObservableType energy {0.0};

  for (unsigned int i = 0; i < Size; i++) {
    if (i != 0) xLeft = i - 1; else xLeft = Size - 1;
    for (unsigned int j = 0; j < Size; j++) {
      if (j != 0) yBelow = j - 1; else yBelow = Size - 1;

        energy += spin[i][j].x * (spin[xLeft][j].x + spin[i][yBelow].x) + 
                  spin[i][j].y * (spin[xLeft][j].y + spin[i][yBelow].y) +
                  spin[i][j].z * (spin[xLeft][j].z + spin[i][yBelow].z);
        
    }  
  }
  
  return -energy; // ferromagnetic (FO) coupling

}


ObservableType Heisenberg2D::getZeemanEnergy()
{
  return 0.0;
}


std::tuple<ObservableType, ObservableType, ObservableType, ObservableType> Heisenberg2D::getMagnetization()
{

  ObservableType m1 {0.0};
  ObservableType m2 {0.0};
  ObservableType m3 {0.0};
  ObservableType m4 {0.0};

  for (unsigned int i = 0; i < Size; i++) {
    for (unsigned int j = 0; j < Size; j++) {
      m1 += spin[i][j].x;
      m2 += spin[i][j].y;
      m3 += spin[i][j].z;
    }
  }

  m4 = sqrt(m1 * m1 + m2 * m2 + m3 * m3);

  return {m1, m2, m3, m4};

}


ObservableType Heisenberg2D::getDifferenceInExchangeInteractions()
{

  unsigned int xLeft, yBelow;
  unsigned int xRight, yAbove;
  ObservableType energyChange {0.0};

  if (CurX != 0) xLeft = CurX - 1; else xLeft = Size - 1;
  if (CurY != 0) yBelow = CurY - 1; else yBelow = Size - 1;
  if (CurX != (Size-1) ) xRight = CurX + 1; else xRight = 0;
  if (CurY != (Size-1) ) yAbove = CurY + 1; else yAbove = 0;

  energyChange = (spin[xLeft][CurY].x + spin[xRight][CurY].x + spin[CurX][yBelow].x + spin[CurX][yAbove].x) * 
                 (spin[CurX][CurY].x - oldSpin.x) +
                 (spin[xLeft][CurY].y + spin[xRight][CurY].y + spin[CurX][yBelow].y + spin[CurX][yAbove].y) * 
                 (spin[CurX][CurY].y - oldSpin.y) +
                 (spin[xLeft][CurY].z + spin[xRight][CurY].z + spin[CurX][yBelow].z + spin[CurX][yAbove].z) * 
                 (spin[CurX][CurY].z - oldSpin.z) ;

  return -energyChange;           // ferromagnetic (FO) coupling

}


ObservableType Heisenberg2D::getDifferenceInZeemanEnergy()
{
  return 0.0;
}


void Heisenberg2D::doMCMove()
{

  double r1, r2, rr;

  // Need this here since resetObservables() is not called if firstTimeGetMeasures = false
  //for (int i = 0; i < numObservables; i++)
  //  oldObservables[i] = observables[i];

  CurX = unsigned(getIntRandomNumber()) % Size;
  CurY = unsigned(getIntRandomNumber()) % Size;

  oldSpin = spin[CurX][CurY];

  do {
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
  spin[CurX][CurY] = oldSpin;
  restoreObservables();
}
*/

void Heisenberg2D::acceptMCMove()
{
  // update "old" observables
  for (unsigned int i=0; i < numObservables; i++)
    oldObservables[i] = observables[i];
}


void Heisenberg2D::rejectMCMove()
{
  spin[CurX][CurY] = oldSpin;
  for (unsigned int i=0; i < numObservables; i++)
    observables[i] = oldObservables[i];
}

/*
void Heisenberg2D::buildMPIConfigurationType()
{
}
*/


void Heisenberg2D::readSpinConfigFile(const std::filesystem::path& spinConfigFile)
{

  std::cout << "\n   Heisenberg2D class reading configuration file: " << spinConfigFile << "\n";

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
            //std::cout << "   Heisenberg2D: numberOfSpins = " << numberOfSpins << "\n";
            continue;
          }
          else if (key == "Observables") {
            unsigned int counter = 0;
            while (lineStream && counter < numObservables) {
              lineStream >> observables[counter];
              //std::cout << "   Heisenberg2D: observables[" << counter << "] = " << observables[counter] << "\n";
              counter++;
            }
            continue;
          }
          else if (key == "SpinConfiguration") {
            //std::cout << "   Heisenberg2D: Spin Configuration read: \n";
            for (unsigned int i=0; i<Size; i++) {
              for (unsigned int j=0; j<Size; j++) {
                lineStream.clear();
                std::getline(inputFile, line);               
                if (!line.empty())  lineStream.str(line);
                lineStream >> spin[i][j].x >> spin[i][j].y >> spin[i][j].z;
                //printf("      %8.5f %8.5f %8.5f\n", spin[i][j][k].x, spin[i][j][k].y, spin[i][j][k].z);
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
  assert(numberOfSpins == systemSize);
  
  printf("   Initial configuration read:\n");
  for (unsigned int i=0; i<Size; i++) {
    for (unsigned int j=0; j<Size; j++)
      printf("      %8.5f %8.5f %8.5f\n", spin[i][j].x, spin[i][j].y, spin[i][j].z);
    printf("\n");
  }

}


void Heisenberg2D::initializeSpinConfiguration(int initial)
{

  double r1, r2, rr;

  for (unsigned int i = 0; i < Size; i++) {
    for (unsigned int j = 0; j < Size; j++) {

      switch (initial) {
        case 1 : {
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
        default  : {
          do {
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
