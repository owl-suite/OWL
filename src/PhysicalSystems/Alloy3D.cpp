#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <filesystem>
#include <fstream>
#include "Alloy3D.hpp"
#include "Utilities/CheckFile.hpp"
#include "Utilities/CompareNumbers.hpp"
#include "Utilities/RandomNumberGenerator.hpp"


// OK
Alloy3D::Alloy3D(const char* inputFile, int initial) : lattice(inputFile)
{

  printf("\n");
  printf("Simulation for alloy with 3D crystal structure: %dx%dx%d unit cells \n", 
         lattice.unitCellDimensions[0], lattice.unitCellDimensions[1], lattice.unitCellDimensions[2]);

  assert (lattice.totalNumberOfAtoms > 0);
  setSystemSize(lattice.totalNumberOfAtoms);
  atom.resize(systemSize);

  if (std::filesystem::exists(inputFile))
    readCompositionInfo(inputFile);

  // Initialize observables
  nearestNeighborPairTypes.resize(numberOfElements, numberOfElements);
  observableName.push_back("Total energy, E");                     // observables[0] : total energy

  char obsName[101];
  for (unsigned int i=0; i<numberOfElements; i++) {
    for (unsigned int j=i; j<numberOfElements; j++) {
      sprintf(obsName, "Prob. pair type %2s-%2s", convertElementToString(elementTypes[i]).c_str(),
                                                  convertElementToString(elementTypes[j]).c_str());
      observableName.push_back(obsName);                           // observables[1-] : nearest neighbor pair types
      //std::cout << i << "  " << j << "  " << obsName << "\n";
    }
  }

  initializeObservables(observableName.size());

  // Initialize configuration from file if applicable
  if (std::filesystem::exists("config_initial.dat"))
    readAtomConfigFile("config_initial.dat");
  else if (simInfo.restartFlag && std::filesystem::exists("configurations/config_checkpoint.dat"))
    readAtomConfigFile("configurations/config_checkpoint.dat");
  else
    initializeAtomConfiguration(initial);

  firstTimeGetMeasures = true;
  getObservablesFromScratch();

  writeConfiguration(0, "configurations/config_initial.dat");

}

// OK
Alloy3D::~Alloy3D()
{

  printf("Exiting Alloy3D class...\n");

}


// OK
void Alloy3D::readCompositionInfo(const std::filesystem::path& mainInputFile)
{

  std::cout << "\n   Alloy3D class reading input file: " << mainInputFile << "\n";

  std::ifstream inputFile(mainInputFile);
  std::string line, key;

  if (inputFile.is_open()) {

    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
        std::istringstream lineStream(line);
        lineStream >> key;

        if (key.compare(0, 1, "#") != 0) {

          if (key == "NumberOfElements") {
            lineStream >> numberOfElements;
            std::cout << "\n     Number of elements in the system = " << numberOfElements << "\n";

            // Allocate memory for other quantities to be read in next
            interactions.resize(numberOfElements, numberOfElements);
            continue;
          }

          else if (key == "ElementTypes") {
            unsigned int counter = 0;
            std::string element;
            while (lineStream && counter < numberOfElements) {
              lineStream >> element;
              elementTypes.push_back(convertStringToElement(element));
              //std::cout << "ElementTypes read = " << element << "\n";
              counter++;
            }
            //std::cout << "ElementTypes counter = " << counter << ", numberOfElements = " << numberOfElements << "\n";
            assert(counter == numberOfElements);

            std::cout << "\n     Element types : ";
            for (auto i : elementTypes)
              std::cout << convertElementToString(i) << "  ";
            std::cout << "\n";
            continue;
          }

          else if (key == "Composition") {
            unsigned int counter {0};
            double comp {0.0};
            while (lineStream && counter < numberOfElements) {
              lineStream >> comp;
              composition.push_back(comp);
              //std::cout << "Composition read = " << comp << "\n";
              counter++;
            }
            assert(counter == numberOfElements);

            std::cout << "\n     Composition : ";
            for (auto i : composition)
              std::cout << i << "  ";
            std::cout << "\n";
            continue;
          }
          
          else if (key == "Interactions") {
            std::cout << "\n     Interactions : ";
            for (unsigned int j=0; j<numberOfElements; j++) {
              lineStream >> interactions(0,j);
              std::cout << interactions(0,j) << "  ";
            }
            std::cout << "\n";

            // Read the rest of the interaction matrix
            for (unsigned int i=1; i<numberOfElements; i++) {
              lineStream.clear();
              std::getline(inputFile, line);               
              if (!line.empty()) lineStream.str(line);
              std::cout << "                    ";
              for (unsigned int j=0; j<numberOfElements; j++) {
                lineStream >> interactions(i,j);
                std::cout << interactions(i,j) << "  ";
              }
              std::cout << "\n";
            }
            continue;
          }

        }
      
      }

    }

    inputFile.close();

  }

}

// OK
// TODO: need to test with a real config file!
void Alloy3D::writeConfiguration(int format, const char* filename)
{

  FILE* configFile;
  if (filename != NULL) configFile = fopen(filename, "w");
  else configFile = stdout;

  switch (format) {

    case 1 : {      // xyz format

      fprintf(configFile, "%u\n", systemSize);
      fprintf(configFile, "FCC ");
      for (unsigned int i=0; i<numberOfElements; i++)
        fprintf(configFile, "X_%2s=%5.3f ", convertElementToString(elementTypes[i]).c_str(), composition[i]);
      fprintf(configFile, "x=%d y=%d z=%d ", lattice.unitCellDimensions[0], lattice.unitCellDimensions[1], lattice.unitCellDimensions[2]);
      fprintf(configFile, "unit_cell=%d \n", systemSize);
      for (unsigned int i = 0; i < systemSize; i++) {
        fprintf(configFile, "%2s %10.6f %10.6f %10.6f\n",
                convertElementToString(atom[i]).c_str(),
                lattice.globalAtomicPositions(0, i), lattice.globalAtomicPositions(1, i), lattice.globalAtomicPositions(2, i));
      }
      break;

    }

    case 0 :
    default : {

      fprintf(configFile, "# Customized 3D alloy structure: %dx%dx%d unit cells \n\n", lattice.unitCellDimensions[0], lattice.unitCellDimensions[1], lattice.unitCellDimensions[2]);
      fprintf(configFile, "TotalNumberOfAtoms %u \n", systemSize);
      fprintf(configFile, "Observables ");
      for (unsigned int i = 0; i < numObservables; i++)
        fprintf(configFile, " %15.8f", observables[i]);
      fprintf(configFile, "\n\n");

      fprintf(configFile, "AlloyConfiguration\n");
      for (unsigned int i = 0; i < systemSize; i++)
        fprintf(configFile, "%s\n", convertElementToString(atom[i]).c_str() );
      fprintf(configFile, "\n");

    }

  }

  if (filename != NULL) fclose(configFile);

}


// OK
void Alloy3D::getObservablesFromScratch()
{

  observables[0] = getExchangeInteractions();
  getAdditionalObservables();
  firstTimeGetMeasures = false;

}


// OK
void Alloy3D::getObservables() 
{

  observables[0] += getDifferenceInExchangeInteractions();
  
}


void Alloy3D::getAdditionalObservables()
{
  
  getNearestNeighborPairTypes();
  unsigned int index = 1;
 
  for (unsigned int i=0; i<numberOfElements; i++) {
    for (unsigned int j=i; j<numberOfElements; j++) {
      if (i==j)
        observables[index] = nearestNeighborPairTypes(i,j);
      else
        observables[index] = 2.0 * nearestNeighborPairTypes(i,j);
      index++;
    }
  }

}


// OK
void Alloy3D::doMCMove()
{

  currentPosition1 = getUnsignedIntRandomNumber() % systemSize;
  currentPosition2 = getUnsignedIntRandomNumber() % systemSize;
  oldAtom1 = atom[currentPosition1];
  oldAtom2 = atom[currentPosition2];

  // Swap the elements of the two chosen sites
  atom[currentPosition1] = oldAtom2;
  atom[currentPosition2] = oldAtom1;


}


// OK
void Alloy3D::acceptMCMove()
{
  // update "old" observables
  for (unsigned int i = 0; i < numObservables; i++)
    oldObservables[i] = observables[i];
}


// OK
void Alloy3D::rejectMCMove()
{
  atom[currentPosition1] = oldAtom1;
  atom[currentPosition2] = oldAtom2;
  for (unsigned int i = 0; i < numObservables; i++)
    observables[i] = oldObservables[i];
}


/*
void Alloy3D::buildMPIConfigurationType()
{
}
*/


// OK
// TODO: need to test with a real config file!
void Alloy3D::readAtomConfigFile(const std::filesystem::path& spinConfigFile)
{
  std::cout << "\n   Alloy3D class reading configuration file: " << spinConfigFile << "\n";

  std::ifstream inputFile(spinConfigFile);
  std::string line, key;
  unsigned int numberOfAtoms {0};

  if (inputFile.is_open()) {

    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
          std::istringstream lineStream(line);
          lineStream >> key;

          if (key.compare(0, 1, "#") != 0) {

            if (key == "TotalNumberOfAtoms") {
              lineStream >> numberOfAtoms;
              assert(numberOfAtoms == systemSize);
              //std::cout << "   Alloy3D: numberOfAtoms = " << numberOfAtoms << "\n";
              continue;
            }
            else if (key == "Observables") {
              unsigned int counter = 0;
              while (lineStream && counter < numObservables) {
                lineStream >> observables[counter];
                //std::cout << "   Alloy3D: observables[" << counter << "] = " << observables[counter] << "\n";
                counter++;
              }
              continue;
            }
            else if (key == "AlloyConfiguration") {
              //std::cout << "   Alloy3D:  Atom configuration read: \n";
              lineStream.clear();
              for (unsigned int atomID=0; atomID<numberOfAtoms; atomID++) {
                std::getline(inputFile, line);               
                //if (!line.empty()) lineStream.str(line);
                //lineStream >> ;
                if (!line.empty()) {
                  lattice.globalAtomicSpecies[atomID] = convertStringToElement(line);
                  atom[atomID] = convertStringToElement(line);
                }
              }
              continue;
            }

          }

      }
    }

    inputFile.close();

  }

  // TODO: check if the read-in configuration agrees with the composition in the input file

}


// OK
// TODO: need to check if it works as expected.
void Alloy3D::initializeAtomConfiguration(int initial)
{

  std::cout << "\n      Initializing atomic configuration...\n";

  // At this point:
  // 1. lattice.globalAtomicSpecies should have been initialized to be the same as lattice.unitCell.atomic_species.
  // 2. readCompositionInfo() should have been called to fill in the composition-related variables/
  // 3. atom is resized to systemSize

  std::vector<unsigned int> currentNumberOfAtoms;     // for each element type
  unsigned int              atomCount {0};
  std::vector<double>       accumulatedComposition;
  bool                      numAtomsOK {false};

  // Determine the number of atoms for each species
  numberOfAtomsForEachElement.resize(numberOfElements);
  currentNumberOfAtoms.resize(numberOfElements);

  accumulatedComposition.resize(numberOfElements);
  accumulatedComposition[0] = composition[0];
  
  for (unsigned int i=0; i<numberOfElements; i++) {
    if (i > 0)
      accumulatedComposition[i] = accumulatedComposition[i-1] + composition[i];
    numberOfAtomsForEachElement[i] = unsigned(round(double(systemSize) * composition[i]));
    atomCount += numberOfAtomsForEachElement[i];
    printf("       - number of %3s : %5u \n", convertElementToString(elementTypes[i]).c_str(), numberOfAtomsForEachElement[i]);
  }
  if (atomCount == systemSize) {
    numAtomsOK = true;
    std::cout << "      Total atom count : " << atomCount << "\n";
  }
  else {
    std::cout << "\n      Fine tuning atomic configuration...\n";
    while (!numAtomsOK) {
      int sign = (systemSize - atomCount) > 0 ? 1 : -1;
      for(unsigned int i=0; i<numberOfElements; i++) {
        numberOfAtomsForEachElement[i] += sign;
        atomCount += sign;
        if (atomCount == systemSize) {
          numAtomsOK = true;
          break;
        }
      }
    }
    for (unsigned int i=0; i<numberOfElements; i++)
      printf("       - number of %3s : %5u \n", convertElementToString(elementTypes[i]).c_str(), numberOfAtomsForEachElement[i]);
    std::cout << "       Adjusted total atom count : " << atomCount << "\n";
  }

  // Either copy element species from unit cell (initial=1), or randomly assign an element type to each atom (initial=0, default)
  switch (initial) {

    case 1 :
    
      std::cout << "\n      Initializing atomic species: copying element type from unit cell... ";

      for (unsigned int atomID=0; atomID<systemSize; atomID++) {
        atom[atomID] = lattice.globalAtomicSpecies[atomID];
        currentNumberOfAtoms[getElementIndex(atom[atomID])]++;
      }

      std::cout << "Done. \n";
      break;
    
    case 0 :
    default :

      std::cout << "\n      Initializing atomic species: randomly assigning an element type to each atom... ";
      unsigned int atomID = 0;
      double r;

      while (atomID < systemSize) {

        r = getRandomNumber2();
        for (unsigned int i=0; i<numberOfElements; i++) {
          if (r < accumulatedComposition[i]) {
            if (currentNumberOfAtoms[i] < numberOfAtomsForEachElement[i]) {
              lattice.globalAtomicSpecies[atomID] = elementTypes[i];
              atom[atomID] = elementTypes[i];
              currentNumberOfAtoms[i]++;
              atomID++;
              break;
            }
            else break;
          }
          else continue;
        }
      }

      std::cout << "Done. \n";
      break;

  }

  // Third, compare and check the configuration with target composition
  for (unsigned int i=0; i<numberOfElements; i++)
    assert(currentNumberOfAtoms[i] == numberOfAtomsForEachElement[i]);

}


// OK
unsigned int Alloy3D::getElementIndex(Element elem){

  unsigned int index {0};

  for (unsigned int i=0; i<numberOfElements; i++) {
    if (elementTypes[i] == elem) {
      index = i;
      break;
    }
  }

  return index;

}


// OK
ObservableType Alloy3D::getExchangeInteractions()
{

  ObservableType energy {0.0};
  unsigned int i, j;

  for (unsigned int atomID=0; atomID<systemSize; atomID++) {
    for (auto neighbor : lattice.neighborList[atomID]) {
      i = getElementIndex(atom[atomID]);
      j = getElementIndex(atom[neighbor.atomID]);
      energy += interactions(i,j);
      //78nearestNeighborPairTypes(i,j)++;
    }
  } 

  return 0.5 * energy;          // the factor of 0.5 is for correcting double counting

}


// OK
// checked that it yields the same energy as from scratch
ObservableType Alloy3D::getDifferenceInExchangeInteractions()
{
  
  ObservableType energyChange  {0.0};
  bool theTwoAtomsAreNeighbors (false);
  unsigned int i, j, k;

  // Check if the two atoms are in each other's neighbor list
  for (auto neighbor : lattice.neighborList[currentPosition1])
    if (neighbor.atomID == currentPosition2)
      theTwoAtomsAreNeighbors = true;
  
  for (auto neighbor : lattice.neighborList[currentPosition2])
    if (neighbor.atomID == currentPosition1)
      theTwoAtomsAreNeighbors = true;

  if (theTwoAtomsAreNeighbors) {           // do not count 
    // First atom
    i = getElementIndex(atom[currentPosition1]);
    j = getElementIndex(oldAtom1);

    for (auto neighbor : lattice.neighborList[currentPosition1]) {
      k = getElementIndex(atom[neighbor.atomID]); 
      if (neighbor.atomID != currentPosition2)
        energyChange += interactions(i,k) - interactions(j,k);
    }

    // Second atom
    i = getElementIndex(atom[currentPosition2]);
    j = getElementIndex(oldAtom2);

    for (auto neighbor : lattice.neighborList[currentPosition2]) {
      k = getElementIndex(atom[neighbor.atomID]); 
      if (neighbor.atomID != currentPosition1)
        energyChange += interactions(i,k) - interactions(j,k);
    }
  }
  else {
    // First atom
    i = getElementIndex(atom[currentPosition1]);
    j = getElementIndex(oldAtom1);

    for (auto neighbor : lattice.neighborList[currentPosition1]) {
      k = getElementIndex(atom[neighbor.atomID]); 
      energyChange += interactions(i,k) - interactions(j,k);
    }

    // Second atom
    i = getElementIndex(atom[currentPosition2]);
    j = getElementIndex(oldAtom2);

    for (auto neighbor : lattice.neighborList[currentPosition2]) {
      k = getElementIndex(atom[neighbor.atomID]); 
      energyChange += interactions(i,k) - interactions(j,k);
    }
  }

  return energyChange;

}



void Alloy3D::getNearestNeighborPairTypes()
{

  unsigned int i, j;
  nearestNeighborPairTypes = 0;

  for (unsigned int atomID=0; atomID<systemSize; atomID++) {
    for (auto neighbor : lattice.neighborList[atomID]) {
      i = getElementIndex(atom[atomID]);
      j = getElementIndex(atom[neighbor.atomID]);
      nearestNeighborPairTypes(i,j)++;
    }
  } 

  // Checks: 
  // 1. the nearestNeighborPairTypes matrix should be symmetrical
  // 2. The matrix elements should add up to (coordination number * number of atoms), which includes double counting
  ObservableType totalPairs {0.0};
  for (i=0; i<numberOfElements; i++) {
    for (j=0; j<numberOfElements; j++) {
      //std::cout << "i : " << i << ", j : " << j << ", nearestNeighborPairTypes(i,j) = " << nearestNeighborPairTypes(i,j) << "\n";
      assert(nearestNeighborPairTypes(i,j) == nearestNeighborPairTypes(j,i));
      totalPairs += nearestNeighborPairTypes(i,j);
    }
  }
  assert(totalPairs == ObservableType(lattice.coordinationNumbers[0] * systemSize));

  // Normalization
  for (i=0; i<numberOfElements; i++) {
    for (j=0; j<numberOfElements; j++)
      nearestNeighborPairTypes(i,j) /= totalPairs;
  }
  
}