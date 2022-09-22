#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include "HeisenbergHexagonal2D.hpp"
#include "Utilities/RandomNumberGenerator.hpp"

HeisenbergHexagonal2D::HeisenbergHexagonal2D(const char* spinConfigFile, int initial)
{

  printf("Simulation for 2D Heisenberg model on a hexagonal lattice: %dx%d \n", simInfo.spinModelLatticeSize, simInfo.spinModelLatticeSize);

  Size = simInfo.spinModelLatticeSize;
  setSystemSize(Size * Size);

  readHamiltonian(simInfo.MCInputFile);
  
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

void HeisenbergHexagonal2D::readHamiltonian(const char* mainInputFile)
{
   //if (GlobalComm.thisMPIrank == 0)
  std::cout << "\n   HeisenbergHexagonal2D class reading input file: " << mainInputFile << "\n\n";

  std::ifstream inputFile(mainInputFile);
  std::string line, key;

  for(int i=0; i<maxShells; i++)
    exchangeParameter[i] = 0.0;
  numShells = 0;
  uniaxialAnisotropy = 0.0;
  externalField[0] = externalField[1] = externalField[2] = 0.0;

  if (inputFile.is_open()) {

    while (std::getline(inputFile, line)) {

	    if (!line.empty()) {
	      std::istringstream lineStream(line);
	      lineStream >> key;

	      if (key.compare(0, 1, "#") != 0) {

		      if (key == "ExchangeInteraction") {
		        while(!lineStream.eof() && numShells < maxShells) {
			        lineStream >> exchangeParameter[numShells];
			        numShells++;
			      }
            continue;
		      }
		      else if (key == "UniaxialAnisotropy") {
		        lineStream >> uniaxialAnisotropy;
            continue;
		      }
		      else if (key == "ExternalField") {
		        lineStream >> externalField[0] >> externalField[1] >> externalField[2];
            continue;
		      }
		    }
	    }
	  }
    
    inputFile.close();
  }

  printf("Exchange Interactions:\n");
  for(int i=0; i<numShells; i++)
    printf("  Shell %d : J = %f\n", i+1, exchangeParameter[i]);
  printf("Uniaxial Anisotropy : K = %f\n", uniaxialAnisotropy);
  printf("External Field : H = (%f, %f, %f)\n\n", externalField[0], externalField[1], externalField[2]);
  
}

HeisenbergHexagonal2D::~HeisenbergHexagonal2D()
{
  for (unsigned int i = 0; i < Size; i++) 
    delete[] spin[i];
  delete[] spin;

  printf("HeisenbergHexagonal2D finished\n");
}


//void HeisenbergHexagonal2D::readCommandLineOptions()
//{ };


void HeisenbergHexagonal2D::writeConfiguration(int format, const char* filename)
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


void HeisenbergHexagonal2D::getObservables()
{

  if (firstTimeGetMeasures) {
    //resetObservables();
    observables[0] = getExchangeInteractions() + getExternalFieldEnergy() + getAnisotropyEnergy();
    std::tie(observables[1], observables[2], observables[3], observables[4]) = getMagnetization();

    firstTimeGetMeasures = false;
    //printf("First time getObservables. \n");
  }
  else {
    observables[0] += getDifferenceInExchangeInteractions() + getDifferenceInExternalFieldEnergy() + getDifferenceInAnisotropyEnergy();
    observables[1] += spin[CurX][CurY].x - CurType.x;
    observables[2] += spin[CurX][CurY].y - CurType.y;
    observables[3] += spin[CurX][CurY].z - CurType.z;
    observables[4] = sqrt(observables[1] * observables[1] + observables[2] * observables[2] + observables[3] * observables[3]);

    //printf("observables = %10.5f %10.5f %10.5f %10.5f %10.5f\n", observables[0], observables[1], observables[2], observables[3], observables[4]);
  }

}

ObservableType HeisenbergHexagonal2D::getShell_1_ExchangeInteractions()
{
  
  unsigned int xPlus, xMinus, yPlus, yMinus;
  ObservableType energy {0.0};
  ObservableType dEnergy;

  // nearest neighbor
  
  for (unsigned int i = 0; i < Size; i++)
    {
    if (i != 0) xMinus = i - 1; else xMinus = Size - 1;
    if (i != Size - 1) xPlus = i + 1; else xPlus = 0;
    for (unsigned int j = 0; j < Size; j++)
      {
	if (j != 0) yMinus = j - 1; else yMinus = Size - 1;
	if (j != Size - 1) yPlus = j + 1; else yPlus = 0;

	// (-1, -1)
	dEnergy = spin[i][j].x * spin[xMinus][yMinus].x +
	  spin[i][j].y * spin[xMinus][yMinus].y +
	  spin[i][j].z * spin[xMinus][yMinus].z;
	// (-1, 0)
	dEnergy += spin[i][j].x * spin[xMinus][j].x +
	  spin[i][j].y * spin[xMinus][j].y +
	  spin[i][j].z * spin[xMinus][j].z;
	// (0, -1)
	dEnergy += spin[i][j].x * spin[i][yMinus].x +
	  spin[i][j].y * spin[i][yMinus].y +
	  spin[i][j].z * spin[i][yMinus].z;
	// (1, 0)
	dEnergy += spin[i][j].x * spin[xPlus][j].x +
	  spin[i][j].y * spin[xPlus][j].y +
	  spin[i][j].z * spin[xPlus][j].z;
	// (0, 1)
	dEnergy += spin[i][j].x * spin[i][yPlus].x +
	  spin[i][j].y * spin[i][yPlus].y +
	  spin[i][j].z * spin[i][yPlus].z;
	// (1, 1)
	dEnergy += spin[i][j].x * spin[xPlus][yPlus].x +
	  spin[i][j].y * spin[xPlus][yPlus].y +
	  spin[i][j].z * spin[xPlus][yPlus].z;
	
        energy += 0.5*exchangeParameter[0]*dEnergy;
    }  
  }
  
  return -energy; // ferromagnetic (FO) coupling
}

ObservableType HeisenbergHexagonal2D::getShell_2_ExchangeInteractions()
{
  
  unsigned int xPlus, xMinus, yPlus, yMinus;
  unsigned int xPlus2, xMinus2, yPlus2, yMinus2;
  ObservableType energy {0.0};
  ObservableType dEnergy;

  // nearest neighbor
  
  for (unsigned int i = 0; i < Size; i++)
    {
      if (i != 0) xMinus = i - 1; else xMinus = Size - 1;
      if (i != Size - 1) xPlus = i + 1; else xPlus = 0;
    
      if (i > 1) xMinus2 = i - 2;
      else if (Size > 1) xMinus2 = Size + i - 2; else xMinus2 = 0;
      if (i < Size - 2) xPlus2 = i + 2;
      else if (Size > 1) xPlus2 = i - Size + 2; else xPlus2 = 0;
    
      for (unsigned int j = 0; j < Size; j++)
	{
	  if (j != 0) yMinus = j - 1; else yMinus = Size - 1;
	  if (j != Size - 1) yPlus = j + 1; else yPlus = 0;
	  
	  if (j > 1) yMinus2 = j - 2;
	  else if (Size > 1) yMinus2 = Size + j - 2; else yMinus2 = 0;
	  if (j < Size - 2) yPlus2 = j + 2;
	  else if (Size > 1) yPlus2 = j - Size + 2; else yPlus2 = 0;

	// (-2, -1)
	dEnergy = spin[i][j].x * spin[xMinus2][yMinus].x +
	  spin[i][j].y * spin[xMinus2][yMinus].y +
	  spin[i][j].z * spin[xMinus2][yMinus].z;
	// (-1, -2)
	dEnergy += spin[i][j].x * spin[xMinus][yMinus2].x +
	  spin[i][j].y * spin[xMinus][yMinus2].y +
	  spin[i][j].z * spin[xMinus][yMinus2].z;
	// (-1, 1)
	dEnergy += spin[i][j].x * spin[xMinus][yPlus].x +
	  spin[i][j].y * spin[xMinus][yPlus].y +
	  spin[i][j].z * spin[xMinus][yPlus].z;
	// (1, -1)
	dEnergy += spin[i][j].x * spin[xPlus][yMinus].x +
	  spin[i][j].y * spin[xPlus][yMinus].y +
	  spin[i][j].z * spin[xPlus][yMinus].z;
	// (1, 2)
	dEnergy += spin[i][j].x * spin[xPlus][yPlus2].x +
	  spin[i][j].y * spin[xPlus][yPlus2].y +
	  spin[i][j].z * spin[xPlus][yPlus2].z;
	// (2, 1)
	dEnergy += spin[i][j].x * spin[xPlus2][yPlus].x +
	  spin[i][j].y * spin[xPlus2][yPlus].y +
	  spin[i][j].z * spin[xPlus2][yPlus].z;
	
        energy += 0.5*exchangeParameter[1]*dEnergy;
    }  
  }
  
  return -energy; // ferromagnetic (FO) coupling
}

ObservableType HeisenbergHexagonal2D::getShell_3_ExchangeInteractions()
{
  
  unsigned int xPlus2, xMinus2, yPlus2, yMinus2;
  ObservableType energy {0.0};
  ObservableType dEnergy;

  // nearest neighbor
  
  for (unsigned int i = 0; i < Size; i++)
    {
      if (i > 1) xMinus2 = i - 2;
      else if (Size > 1) xMinus2 = Size + i - 2; else xMinus2 = 0;
      if (i < Size - 2) xPlus2 = i + 2;
      else if (Size > 1) xPlus2 = i - Size + 2; else xPlus2 = 0;
      
      for (unsigned int j = 0; j < Size; j++)
	{
	  if (j > 1) yMinus2 = j - 2;
	  else if (Size > 1) yMinus2 = Size + j - 2; else yMinus2 = 0;
	  if (j < Size - 2) yPlus2 = j + 2;
	  else if (Size > 1) yPlus2 = j - Size + 2; else yPlus2 = 0;

	  // (-2, -2)
	  dEnergy = spin[i][j].x * spin[xMinus2][yMinus2].x +
	    spin[i][j].y * spin[xMinus2][yMinus2].y +
	    spin[i][j].z * spin[xMinus2][yMinus2].z;
	  // (-2, 0)
	  dEnergy += spin[i][j].x * spin[xMinus2][j].x +
	    spin[i][j].y * spin[xMinus2][j].y +
	    spin[i][j].z * spin[xMinus2][j].z;
	  // (0, -2)
	  dEnergy += spin[i][j].x * spin[i][yMinus2].x +
	    spin[i][j].y * spin[i][yMinus2].y +
	    spin[i][j].z * spin[i][yMinus2].z;
	  // (2, 0)
	  dEnergy += spin[i][j].x * spin[xPlus2][j].x +
	    spin[i][j].y * spin[xPlus2][j].y +
	    spin[i][j].z * spin[xPlus2][j].z;
	  // (0, 2)
	  dEnergy += spin[i][j].x * spin[i][yPlus2].x +
	    spin[i][j].y * spin[i][yPlus2].y +
	    spin[i][j].z * spin[i][yPlus2].z;
	  // (2, 2)
	  dEnergy += spin[i][j].x * spin[xPlus2][yPlus2].x +
	    spin[i][j].y * spin[xPlus2][yPlus2].y +
	    spin[i][j].z * spin[xPlus2][yPlus2].z;
	  
	  energy += 0.5*exchangeParameter[2]*dEnergy;
	}  
    }
  
  return -energy; // ferromagnetic (FO) coupling
}

ObservableType HeisenbergHexagonal2D::getShell_4_ExchangeInteractions()
{
  
  unsigned int xPlus, xMinus, yPlus, yMinus;
  unsigned int xPlus2, xMinus2, yPlus2, yMinus2;
  unsigned int xPlus3, xMinus3, yPlus3, yMinus3;
  ObservableType energy {0.0};
  ObservableType dEnergy;

  // nearest neighbor
  
  for (unsigned int i = 0; i < Size; i++)
    {
      if (i != 0) xMinus = i - 1; else xMinus = Size - 1;
      if (i != Size - 1) xPlus = i + 1; else xPlus = 0;
    
      if (i > 1) xMinus2 = i - 2;
      else if (Size > 1) xMinus2 = Size + i - 2; else xMinus2 = 0;
      if (i < Size - 2) xPlus2 = i + 2;
      else if (Size > 1) xPlus2 = i - Size + 2; else xPlus2 = 0;

      if (i > 2) xMinus3 = i - 3;
      else if (Size > 2) xMinus3 = Size - 3 + i;
      else if (Size == 2) xMinus3 = 1; else xMinus3 = 0;    
      if (i < Size - 3) xPlus3 = i + 3;
      else if (Size > 2) xPlus3 = i - Size + 3;
      else if (Size > 1) xPlus3 = 1; else xPlus3 = 0;
    
      for (unsigned int j = 0; j < Size; j++)
	{
	  if (j != 0) yMinus = j - 1; else yMinus = Size - 1;
	  if (j != Size - 1) yPlus = j + 1; else yPlus = 0;
	  
	  if (j > 1) yMinus2 = j - 2;
	  else if (Size > 1) yMinus2 = Size + j - 2; else yMinus2 = 0;
	  if (j < Size - 2) yPlus2 = j + 2;
	  else if (Size > 1) yPlus2 = j - Size + 2; else yPlus2 = 0;

	  if (j > 2) yMinus3 = j - 3;
	  else if (Size > 2) yMinus3 = Size - 3 + j;
	  else if (Size == 2) yMinus3 = 1; else yMinus3 = 0;    
	  if (j < Size - 3) yPlus3 = j + 3;
	  else if (Size > 2) yPlus3 = j - Size + 3;
	  else if (Size > 1) yPlus3 = 1; else yPlus3 = 0;

	// (-3, -2)
	dEnergy = spin[i][j].x * spin[xMinus3][yMinus2].x +
	  spin[i][j].y * spin[xMinus3][yMinus2].y +
	  spin[i][j].z * spin[xMinus3][yMinus2].z;
	// (-3, -1)
	dEnergy += spin[i][j].x * spin[xMinus3][yMinus].x +
	  spin[i][j].y * spin[xMinus3][yMinus].y +
	  spin[i][j].z * spin[xMinus3][yMinus].z;
	// (-2, -3)
	dEnergy += spin[i][j].x * spin[xMinus2][yMinus3].x +
	  spin[i][j].y * spin[xMinus2][yMinus3].y +
	  spin[i][j].z * spin[xMinus2][yMinus3].z;
	// (-2, 1)
	dEnergy += spin[i][j].x * spin[xMinus2][yPlus].x +
	  spin[i][j].y * spin[xMinus2][yPlus].y +
	  spin[i][j].z * spin[xMinus2][yPlus].z;
	// (-1, -3)
	dEnergy += spin[i][j].x * spin[xMinus][yMinus3].x +
	  spin[i][j].y * spin[xMinus][yMinus3].y +
	  spin[i][j].z * spin[xMinus][yMinus3].z;
	// (-1, 2)
	dEnergy += spin[i][j].x * spin[xMinus][yPlus2].x +
	  spin[i][j].y * spin[xMinus][yPlus2].y +
	  spin[i][j].z * spin[xMinus][yPlus2].z;
	// (1, -2)
	dEnergy = spin[i][j].x * spin[xPlus][yMinus2].x +
	  spin[i][j].y * spin[xPlus][yMinus2].y +
	  spin[i][j].z * spin[xPlus][yMinus2].z;
	// (1, 3)
	dEnergy += spin[i][j].x * spin[xPlus][yPlus3].x +
	  spin[i][j].y * spin[xPlus][yPlus3].y +
	  spin[i][j].z * spin[xPlus][yPlus3].z;
	// (2, -1)
	dEnergy += spin[i][j].x * spin[xPlus2][yMinus].x +
	  spin[i][j].y * spin[xPlus2][yMinus].y +
	  spin[i][j].z * spin[xPlus2][yMinus].z;
	// (2, 3)
	dEnergy += spin[i][j].x * spin[xPlus2][yPlus3].x +
	  spin[i][j].y * spin[xPlus2][yPlus3].y +
	  spin[i][j].z * spin[xPlus2][yPlus3].z;
	// (3, 1)
	dEnergy += spin[i][j].x * spin[xPlus3][yPlus].x +
	  spin[i][j].y * spin[xPlus3][yPlus].y +
	  spin[i][j].z * spin[xPlus3][yPlus].z;
	// (3, 2)
	dEnergy += spin[i][j].x * spin[xPlus3][yPlus2].x +
	  spin[i][j].y * spin[xPlus3][yPlus2].y +
	  spin[i][j].z * spin[xPlus3][yPlus2].z;
	
        energy += 0.5*exchangeParameter[3]*dEnergy;
    }  
  }
  
  return -energy; // ferromagnetic (FO) coupling
}

ObservableType HeisenbergHexagonal2D::getShell_5_ExchangeInteractions()
{
  
  unsigned int xPlus3, xMinus3, yPlus3, yMinus3;
  ObservableType energy {0.0};
  ObservableType dEnergy;

  // nearest neighbor
  
  for (unsigned int i = 0; i < Size; i++)
    {
      if (i > 2) xMinus3 = i - 3;
      else if (Size > 2) xMinus3 = Size - 3 + i;
      else if (Size == 2) xMinus3 = 1; else xMinus3 = 0;    
      if (i < Size - 3) xPlus3 = i + 3;
      else if (Size > 2) xPlus3 = i - Size + 3;
      else if (Size > 1) xPlus3 = 1; else xPlus3 = 0;
      
      for (unsigned int j = 0; j < Size; j++)
      {
	if (j > 2) yMinus3 = j - 3;
	else if (Size > 2) yMinus3 = Size - 3 + j;
	else if (Size == 2) yMinus3 = 1; else yMinus3 = 0;    
	if (j < Size - 3) yPlus3 = j + 3;
	else if (Size > 2) yPlus3 = j - Size + 3;
	else if (Size > 1) yPlus3 = 1; else yPlus3 = 0;

	// (-3, -3)
	dEnergy = spin[i][j].x * spin[xMinus3][yMinus3].x +
	  spin[i][j].y * spin[xMinus3][yMinus3].y +
	  spin[i][j].z * spin[xMinus3][yMinus3].z;
	// (-3, 0)
	dEnergy += spin[i][j].x * spin[xMinus3][j].x +
	  spin[i][j].y * spin[xMinus3][j].y +
	  spin[i][j].z * spin[xMinus3][j].z;
	// (0, -3)
	dEnergy += spin[i][j].x * spin[i][yMinus3].x +
	  spin[i][j].y * spin[i][yMinus3].y +
	  spin[i][j].z * spin[i][yMinus3].z;
	// (3, 0)
	dEnergy += spin[i][j].x * spin[xPlus3][j].x +
	  spin[i][j].y * spin[xPlus3][j].y +
	  spin[i][j].z * spin[xPlus3][j].z;
	// (0, 3)
	dEnergy += spin[i][j].x * spin[i][yPlus3].x +
	  spin[i][j].y * spin[i][yPlus3].y +
	  spin[i][j].z * spin[i][yPlus3].z;
	// (3, 3)
	dEnergy += spin[i][j].x * spin[xPlus3][yPlus3].x +
	  spin[i][j].y * spin[xPlus3][yPlus3].y +
	  spin[i][j].z * spin[xPlus3][yPlus3].z;
	
        energy += 0.5*exchangeParameter[4]*dEnergy;
    }  
  }
  
  return -energy; // ferromagnetic (FO) coupling
}

ObservableType HeisenbergHexagonal2D::getExchangeInteractions()
{

  ObservableType energy {0.0};
  
  switch(numShells)
    {
    case 5: energy += getShell_5_ExchangeInteractions();
    case 4: energy += getShell_4_ExchangeInteractions();
    case 3: energy += getShell_3_ExchangeInteractions();
    case 2: energy += getShell_2_ExchangeInteractions();
    case 1: energy += getShell_1_ExchangeInteractions();
    }
  
  return energy;

}


ObservableType HeisenbergHexagonal2D::getExternalFieldEnergy()
{
  ObservableType energy {0.0};

    for (unsigned int i = 0; i < Size; i++)
      for (unsigned int j = 0; j < Size; j++)
        energy += spin[i][j].x * externalField[0] + 
                  spin[i][j].y * externalField[1] +
                  spin[i][j].z * externalField[2];

  return -energy;
}

ObservableType HeisenbergHexagonal2D::getAnisotropyEnergy()
{
  ObservableType energy {0.0};

    for (unsigned int i = 0; i < Size; i++)
      for (unsigned int j = 0; j < Size; j++)
        energy += spin[i][j].z * spin[i][j].z; 

  return uniaxialAnisotropy * energy;
}


std::tuple<ObservableType, ObservableType, ObservableType, ObservableType> HeisenbergHexagonal2D::getMagnetization()
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

ObservableType HeisenbergHexagonal2D::getShell_1_DifferenceInExchangeInteractions()
{
  
  unsigned int xPlus, xMinus, yPlus, yMinus;

  ObservableType dEnergy;

  SpinDirection dCur;
  
  int i = CurX;
  int j = CurY;

  dCur.x = spin[CurX][CurY].x - CurType.x;
  dCur.y = spin[CurX][CurY].y - CurType.y;
  dCur.z = spin[CurX][CurY].z - CurType.z;

  if (i != 0) xMinus = i - 1; else xMinus = Size - 1;
  if (i != Size - 1) xPlus = i + 1; else xPlus = 0;

  if (j != 0) yMinus = j - 1; else yMinus = Size - 1;
  if (j != Size - 1) yPlus = j + 1; else yPlus = 0;

  // (-1, -1)
  dEnergy = dCur.x * spin[xMinus][yMinus].x +
    dCur.y * spin[xMinus][yMinus].y +
    dCur.z * spin[xMinus][yMinus].z;
  // (-1, 0)
  dEnergy += dCur.x * spin[xMinus][j].x +
    dCur.y * spin[xMinus][j].y +
    dCur.z * spin[xMinus][j].z;
  // (0, -1)
  dEnergy += dCur.x * spin[i][yMinus].x +
    dCur.y * spin[i][yMinus].y +
    dCur.z * spin[i][yMinus].z;
  // (1, 0)
  dEnergy += dCur.x * spin[xPlus][j].x +
    dCur.y * spin[xPlus][j].y +
    dCur.z * spin[xPlus][j].z;
  // (0, 1)
  dEnergy += dCur.x * spin[i][yPlus].x +
    dCur.y * spin[i][yPlus].y +
    dCur.z * spin[i][yPlus].z;
  // (1, 1)
  dEnergy += dCur.x * spin[xPlus][yPlus].x +
    dCur.y * spin[xPlus][yPlus].y +
    dCur.z * spin[xPlus][yPlus].z;
	
  //     energy += 0.5*exchangeParameter[0]*dEnergy;
  
  return -exchangeParameter[0]*dEnergy; // ferromagnetic (FO) coupling
}

ObservableType HeisenbergHexagonal2D::getShell_2_DifferenceInExchangeInteractions()
{
  
  unsigned int xPlus, xMinus, yPlus, yMinus;
  unsigned int xPlus2, xMinus2, yPlus2, yMinus2;

  ObservableType dEnergy;

  SpinDirection dCur;
  
  int i = CurX;
  int j = CurY;

  dCur.x = spin[CurX][CurY].x - CurType.x;
  dCur.y = spin[CurX][CurY].y - CurType.y;
  dCur.z = spin[CurX][CurY].z - CurType.z;
  

  if (i != 0) xMinus = i - 1; else xMinus = Size - 1;
  if (i != Size - 1) xPlus = i + 1; else xPlus = 0;
    
  if (i > 1) xMinus2 = i - 2;
  else if (Size > 1) xMinus2 = Size + i - 2; else xMinus2 = 0;
  if (i < Size - 2) xPlus2 = i + 2;
  else if (Size > 1) xPlus2 = i - Size + 2; else xPlus2 = 0;
    
  if (j != 0) yMinus = j - 1; else yMinus = Size - 1;
  if (j != Size - 1) yPlus = j + 1; else yPlus = 0;
	  
  if (j > 1) yMinus2 = j - 2;
  else if (Size > 1) yMinus2 = Size + j - 2; else yMinus2 = 0;
  if (j < Size - 2) yPlus2 = j + 2;
  else if (Size > 1) yPlus2 = j - Size + 2; else yPlus2 = 0;

  // (-2, -1)
  dEnergy = dCur.x * spin[xMinus2][yMinus].x +
    dCur.y * spin[xMinus2][yMinus].y +
    dCur.z * spin[xMinus2][yMinus].z;
  // (-1, -2)
  dEnergy += dCur.x * spin[xMinus][yMinus2].x +
    dCur.y * spin[xMinus][yMinus2].y +
    dCur.z * spin[xMinus][yMinus2].z;
  // (-1, 1)
  dEnergy += dCur.x * spin[xMinus][yPlus].x +
    dCur.y * spin[xMinus][yPlus].y +
    dCur.z * spin[xMinus][yPlus].z;
  // (1, -1)
  dEnergy += dCur.x * spin[xPlus][yMinus].x +
    dCur.y * spin[xPlus][yMinus].y +
    dCur.z * spin[xPlus][yMinus].z;
  // (1, 2)
  dEnergy += dCur.x * spin[xPlus][yPlus2].x +
    dCur.y * spin[xPlus][yPlus2].y +
    dCur.z * spin[xPlus][yPlus2].z;
  // (2, 1)
  dEnergy += dCur.x * spin[xPlus2][yPlus].x +
    dCur.y * spin[xPlus2][yPlus].y +
    dCur.z * spin[xPlus2][yPlus].z;
  
  // energy += 0.5*exchangeParameter[1]*dEnergy;
  
  return -exchangeParameter[1]*dEnergy; // ferromagnetic (FO) coupling
}

ObservableType HeisenbergHexagonal2D::getShell_3_DifferenceInExchangeInteractions()
{
  
  unsigned int xPlus2, xMinus2, yPlus2, yMinus2;

  ObservableType dEnergy;

  SpinDirection dCur;
  
  int i = CurX;
  int j = CurY;

  dCur.x = spin[CurX][CurY].x - CurType.x;
  dCur.y = spin[CurX][CurY].y - CurType.y;
  dCur.z = spin[CurX][CurY].z - CurType.z;
  
  if (i > 1) xMinus2 = i - 2;
  else if (Size > 1) xMinus2 = Size + i - 2; else xMinus2 = 0;
  if (i < Size - 2) xPlus2 = i + 2;
  else if (Size > 1) xPlus2 = i - Size + 2; else xPlus2 = 0;
      
  if (j > 1) yMinus2 = j - 2;
  else if (Size > 1) yMinus2 = Size + j - 2; else yMinus2 = 0;
  if (j < Size - 2) yPlus2 = j + 2;
  else if (Size > 1) yPlus2 = j - Size + 2; else yPlus2 = 0;

  // (-2, -2)
  dEnergy = dCur.x * spin[xMinus2][yMinus2].x +
    dCur.y * spin[xMinus2][yMinus2].y +
    dCur.z * spin[xMinus2][yMinus2].z;
  // (-2, 0)
  dEnergy += dCur.x * spin[xMinus2][j].x +
    dCur.y * spin[xMinus2][j].y +
    dCur.z * spin[xMinus2][j].z;
  // (0, -2)
  dEnergy += dCur.x * spin[i][yMinus2].x +
    dCur.y * spin[i][yMinus2].y +
    dCur.z * spin[i][yMinus2].z;
  // (2, 0)
  dEnergy += dCur.x * spin[xPlus2][j].x +
    dCur.y * spin[xPlus2][j].y +
    dCur.z * spin[xPlus2][j].z;
  // (0, 2)
  dEnergy += dCur.x * spin[i][yPlus2].x +
    dCur.y * spin[i][yPlus2].y +
    dCur.z * spin[i][yPlus2].z;
  // (2, 2)
  dEnergy += dCur.x * spin[xPlus2][yPlus2].x +
    dCur.y * spin[xPlus2][yPlus2].y +
    dCur.z * spin[xPlus2][yPlus2].z;
  
  // energy += 0.5*exchangeParameter[2]*dEnergy;
  
  return -exchangeParameter[2]*dEnergy; // ferromagnetic (FO) coupling
}

ObservableType HeisenbergHexagonal2D::getShell_4_DifferenceInExchangeInteractions()
{
  
  unsigned int xPlus, xMinus, yPlus, yMinus;
  unsigned int xPlus2, xMinus2, yPlus2, yMinus2;
  unsigned int xPlus3, xMinus3, yPlus3, yMinus3;

  ObservableType dEnergy;

  SpinDirection dCur;
  
  int i = CurX;
  int j = CurY;

  dCur.x = spin[CurX][CurY].x - CurType.x;
  dCur.y = spin[CurX][CurY].y - CurType.y;
  dCur.z = spin[CurX][CurY].z - CurType.z;
  
  if (i != 0) xMinus = i - 1; else xMinus = Size - 1;
  if (i != Size - 1) xPlus = i + 1; else xPlus = 0;
    
  if (i > 1) xMinus2 = i - 2;
  else if (Size > 1) xMinus2 = Size + i - 2; else xMinus2 = 0;
  if (i < Size - 2) xPlus2 = i + 2;
  else if (Size > 1) xPlus2 = i - Size + 2; else xPlus2 = 0;

  if (i > 2) xMinus3 = i - 3;
  else if (Size > 2) xMinus3 = Size - 3 + i;
  else if (Size == 2) xMinus3 = 1; else xMinus3 = 0;    
  if (i < Size - 3) xPlus3 = i + 3;
  else if (Size > 2) xPlus3 = i - Size + 3;
  else if (Size > 1) xPlus3 = 1; else xPlus3 = 0;
    
  if (j != 0) yMinus = j - 1; else yMinus = Size - 1;
  if (j != Size - 1) yPlus = j + 1; else yPlus = 0;
	  
  if (j > 1) yMinus2 = j - 2;
  else if (Size > 1) yMinus2 = Size + j - 2; else yMinus2 = 0;
  if (j < Size - 2) yPlus2 = j + 2;
  else if (Size > 1) yPlus2 = j - Size + 2; else yPlus2 = 0;

  if (j > 2) yMinus3 = j - 3;
  else if (Size > 2) yMinus3 = Size - 3 + j;
  else if (Size == 2) yMinus3 = 1; else yMinus3 = 0;    
  if (j < Size - 3) yPlus3 = j + 3;
  else if (Size > 2) yPlus3 = j - Size + 3;
  else if (Size > 1) yPlus3 = 1; else yPlus3 = 0;

  // (-3, -2)
  dEnergy = dCur.x * spin[xMinus3][yMinus2].x +
    dCur.y * spin[xMinus3][yMinus2].y +
    dCur.z * spin[xMinus3][yMinus2].z;
  // (-3, -1)
  dEnergy += dCur.x * spin[xMinus3][yMinus].x +
    dCur.y * spin[xMinus3][yMinus].y +
    dCur.z * spin[xMinus3][yMinus].z;
  // (-2, -3)
  dEnergy += dCur.x * spin[xMinus2][yMinus3].x +
    dCur.y * spin[xMinus2][yMinus3].y +
    dCur.z * spin[xMinus2][yMinus3].z;
  // (-2, 1)
  dEnergy += dCur.x * spin[xMinus2][yPlus].x +
    dCur.y * spin[xMinus2][yPlus].y +
    dCur.z * spin[xMinus2][yPlus].z;
  // (-1, -3)
  dEnergy += dCur.x * spin[xMinus][yMinus3].x +
    dCur.y * spin[xMinus][yMinus3].y +
    dCur.z * spin[xMinus][yMinus3].z;
  // (-1, 2)
  dEnergy += dCur.x * spin[xMinus][yPlus2].x +
    dCur.y * spin[xMinus][yPlus2].y +
    dCur.z * spin[xMinus][yPlus2].z;
  // (1, -2)
  dEnergy = dCur.x * spin[xPlus][yMinus2].x +
    dCur.y * spin[xPlus][yMinus2].y +
    dCur.z * spin[xPlus][yMinus2].z;
  // (1, 3)
  dEnergy += dCur.x * spin[xPlus][yPlus3].x +
    dCur.y * spin[xPlus][yPlus3].y +
    dCur.z * spin[xPlus][yPlus3].z;
  // (2, -1)
  dEnergy += dCur.x * spin[xPlus2][yMinus].x +
    dCur.y * spin[xPlus2][yMinus].y +
    dCur.z * spin[xPlus2][yMinus].z;
  // (2, 3)
  dEnergy += dCur.x * spin[xPlus2][yPlus3].x +
    dCur.y * spin[xPlus2][yPlus3].y +
    dCur.z * spin[xPlus2][yPlus3].z;
  // (3, 1)
  dEnergy += dCur.x * spin[xPlus3][yPlus].x +
    dCur.y * spin[xPlus3][yPlus].y +
    dCur.z * spin[xPlus3][yPlus].z;
  // (3, 2)
  dEnergy += dCur.x * spin[xPlus3][yPlus2].x +
    dCur.y * spin[xPlus3][yPlus2].y +
    dCur.z * spin[xPlus3][yPlus2].z;
  
  // energy += 0.5*exchangeParameter[3]*dEnergy;
  
  return -exchangeParameter[3]*dEnergy; // ferromagnetic (FO) coupling
}

ObservableType HeisenbergHexagonal2D::getShell_5_DifferenceInExchangeInteractions()
{
  
  unsigned int xPlus3, xMinus3, yPlus3, yMinus3;

  ObservableType dEnergy;

  SpinDirection dCur;
  
  int i = CurX;
  int j = CurY;

  dCur.x = spin[CurX][CurY].x - CurType.x;
  dCur.y = spin[CurX][CurY].y - CurType.y;
  dCur.z = spin[CurX][CurY].z - CurType.z;
  
  if (i > 2) xMinus3 = i - 3;
  else if (Size > 2) xMinus3 = Size - 3 + i;
  else if (Size == 2) xMinus3 = 1; else xMinus3 = 0;    
  if (i < Size - 3) xPlus3 = i + 3;
  else if (Size > 2) xPlus3 = i - Size + 3;
  else if (Size > 1) xPlus3 = 1; else xPlus3 = 0;
      

  if (j > 2) yMinus3 = j - 3;
  else if (Size > 2) yMinus3 = Size - 3 + j;
  else if (Size == 2) yMinus3 = 1; else yMinus3 = 0;    
  if (j < Size - 3) yPlus3 = j + 3;
  else if (Size > 2) yPlus3 = j - Size + 3;
  else if (Size > 1) yPlus3 = 1; else yPlus3 = 0;

  // (-3, -3)
  dEnergy = dCur.x * spin[xMinus3][yMinus3].x +
    dCur.y * spin[xMinus3][yMinus3].y +
    dCur.z * spin[xMinus3][yMinus3].z;
  // (-3, 0)
  dEnergy += dCur.x * spin[xMinus3][j].x +
    dCur.y * spin[xMinus3][j].y +
    dCur.z * spin[xMinus3][j].z;
  // (0, -3)
  dEnergy += dCur.x * spin[i][yMinus3].x +
    dCur.y * spin[i][yMinus3].y +
    dCur.z * spin[i][yMinus3].z;
  // (3, 0)
  dEnergy += dCur.x * spin[xPlus3][j].x +
    dCur.y * spin[xPlus3][j].y +
    dCur.z * spin[xPlus3][j].z;
  // (0, 3)
  dEnergy += dCur.x * spin[i][yPlus3].x +
    dCur.y * spin[i][yPlus3].y +
    dCur.z * spin[i][yPlus3].z;
  // (3, 3)
  dEnergy += dCur.x * spin[xPlus3][yPlus3].x +
    dCur.y * spin[xPlus3][yPlus3].y +
    dCur.z * spin[xPlus3][yPlus3].z;
	
  // energy += 0.5*exchangeParameter[4]*dEnergy;

  return -exchangeParameter[4]*dEnergy; // ferromagnetic (FO) coupling
}

ObservableType HeisenbergHexagonal2D::getDifferenceInExchangeInteractions()
{

  ObservableType energyChange {0.0};
  
  switch(numShells)
    {
    case 5: energyChange += getShell_5_DifferenceInExchangeInteractions();
    case 4: energyChange += getShell_4_DifferenceInExchangeInteractions();
    case 3: energyChange += getShell_3_DifferenceInExchangeInteractions();
    case 2: energyChange += getShell_2_DifferenceInExchangeInteractions();
    case 1: energyChange += getShell_1_DifferenceInExchangeInteractions();
    }
  
  return energyChange;

}


ObservableType HeisenbergHexagonal2D::getDifferenceInExternalFieldEnergy()
{
  ObservableType energyChange {0.0};

  energyChange = (spin[CurX][CurY].x - CurType.x) * externalField[0] +
                 (spin[CurX][CurY].y - CurType.y) * externalField[1] +
                 (spin[CurX][CurY].z - CurType.z) * externalField[2];
  
  return energyChange;
}

ObservableType HeisenbergHexagonal2D::getDifferenceInAnisotropyEnergy()
{
  ObservableType energyChange {0.0};

  energyChange = spin[CurX][CurY].z * spin[CurX][CurY].z - CurType.z * CurType.z;
  
  return uniaxialAnisotropy * energyChange;
}

void HeisenbergHexagonal2D::doMCMove()
{

  double r1, r2, rr;

  // Need this here since resetObservables() is not called if firstTimeGetMeasures = false
  //for (int i = 0; i < numObservables; i++)
  //  oldObservables[i] = observables[i];

  CurX = unsigned(getIntRandomNumber()) % Size;
  CurY = unsigned(getIntRandomNumber()) % Size;

  CurType = spin[CurX][CurY];

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
void HeisenbergHexagonal2D::undoMCMove()
{
  spin[CurX][CurY] = CurType;
  restoreObservables();
}
*/

void HeisenbergHexagonal2D::acceptMCMove()
{
  // update "old" observables
  for (unsigned int i=0; i < numObservables; i++)
    oldObservables[i] = observables[i];
}


void HeisenbergHexagonal2D::rejectMCMove()
{
  spin[CurX][CurY] = CurType;
  for (unsigned int i=0; i < numObservables; i++)
    observables[i] = oldObservables[i];
}

/*
void HeisenbergHexagonal2D::buildMPIConfigurationType()
{
}
*/


void HeisenbergHexagonal2D::readSpinConfigFile(const std::filesystem::path& spinConfigFile)
{

  std::cout << "\n   HeisenbergHexagonal2D class reading configuration file: " << spinConfigFile << "\n";

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
            //std::cout << "   HeisenbergHexagonal2D: numberOfSpins = " << numberOfSpins << "\n";
            continue;
          }
          else if (key == "Observables") {
            unsigned int counter = 0;
            while (lineStream && counter < numObservables) {
              lineStream >> observables[counter];
              //std::cout << "   HeisenbergHexagonal2D: observables[" << counter << "] = " << observables[counter] << "\n";
              counter++;
            }
            continue;
          }
          else if (key == "SpinConfiguration") {
            //std::cout << "   HeisenbergHexagonal2D: Spin Configuration read: \n";
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


void HeisenbergHexagonal2D::initializeSpinConfiguration(int initial)
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
