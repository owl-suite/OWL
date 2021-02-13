#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <vector>
#include "CrystalBase.hpp"
#include "Utilities/CompareNumbers.hpp"


// Note (data layout of Matrix class):
// (i) mat(i,j)
//     resize(nRow,nCol); nRow = dim. of coordinates, nCol = # of atoms
//       j ->
//     i    x0  x1 ... x[j-1]
//     |    y0  y1     y[j-1]
//     v    z0  z1 ... z[j-1]
//
// (ii) mat[k]
//      index is calculated by k = j*lDim + i
//       0  1  2  3  4  5 ... 
//      x0 y0 z0 x1 y1 z1 ... 


// Constructor 1: initialize unit cell and lattice from input file
Lattice::Lattice(const char* inputFile)
{
  
  std::cout << "\nInitializing crystal structure... \n";

  // Initialize unit cell
  unitCell.lattice_vectors.resize(3, 3);
  if (std::filesystem::exists(inputFile))
    readUnitCellInfo(inputFile);
  else {
    std::cout << "Input file " << inputFile << " does not exist! \n";
    std::cout << "OWL aborting...\n";
    exit(7);
  }
  
  numberOfUnitCells = unitCellDimensions[0] * unitCellDimensions[1] * unitCellDimensions[2];
  constructUnitCellVectors();

  // Initialize crystal
  totalNumberOfAtoms = numberOfUnitCells * unitCell.number_of_atoms;
  constructGlobalCoordinates();
  initializeAtomicSpecies();

  writeAtomicPositions();
  writeAtomicPositions("configurations/atomic_positions.dat");
  // Check:
  //printAllPairwiseDistances();


  // Initialize nearest neighbor lists for each unit cell
  numAdjacentUnitCells = 1;                                   // TODO: to be read in from input file. Default should be all unit cells
  constructRelativeUnitCellVectors();
  nearestNeighborUnitCellList.resize(numberOfUnitCells);
  for (unsigned int i=0; i<numberOfUnitCells; i++)
    nearestNeighborUnitCellList[i] = constructNearestNeighborUnitCellList(i);

  constructRelativeCoordinates();

  // Check:
  //printPairwiseDistancesInUnitCellList(13);

  // Initialize nearest neighbor lists for each atom
  neighborList.resize(totalNumberOfAtoms);
  constructPrimaryNeighborList();
  mapPrimaryToAllNeighborLists();
  
  //for (unsigned int i=0; i<systemSize; i++)
  //  neighborList[i] = constructNeighborListFromNeighboringUnitCells(i);

}


// TODO: Constructor for initializing unit cell in derived classes where crystal structures are known and standard 



// Member functions:

// Note: unit cell layout in lattice_vectors(i,j):
//       j ->
//     i    ax  bx  cx
//     |    ay  by  cy
//     v    az  bz  cz
void Lattice::readUnitCellInfo(const char* mainInputFile)
{

    //if (GlobalComm.thisMPIrank == 0)
      std::cout << "\n   Lattice class reading input file: " << mainInputFile << "\n\n";

    std::ifstream inputFile(mainInputFile);
    std::string line, key;
    unsigned int atom_counter {0};

    if (inputFile.is_open()) {

      while (std::getline(inputFile, line)) {

        //std::cout << "Debug: " << line << "\n"; 

        if (!line.empty()) {
          std::istringstream lineStream(line);
          lineStream >> key;

          //std::cout << "Debug: " << key << "\n"; 

          if (key.compare(0, 1, "#") != 0) {

            if (key == "NumberOfAtomsPerUnitCell") {
              lineStream >> unitCell.number_of_atoms;
              std::cout << "\n     Number of atoms per unit cell = " << unitCell.number_of_atoms << "\n";

              // Allocate memory for atomic_positions
              unitCell.atomic_positions.resize(3, unitCell.number_of_atoms);
              continue;
            }
            else if (key == "LatticeVectors") {
              lineStream >> unitCell.lattice_vectors(0,0) >> unitCell.lattice_vectors(1,0) >> unitCell.lattice_vectors(2,0);
              std::cout << "     Lattice vectors: (" 
                        << unitCell.lattice_vectors(0,0) << ", " 
                        << unitCell.lattice_vectors(1,0) << ", " 
                        << unitCell.lattice_vectors(2,0) << ") \n";

              // Read the other two lattice vectors
              for (unsigned int i=1; i<3; i++) {
                lineStream.clear();
                std::getline(inputFile, line);               
                if (!line.empty()) lineStream.str(line);
                lineStream >> unitCell.lattice_vectors(0,i) >> unitCell.lattice_vectors(1,i) >> unitCell.lattice_vectors(2,i);
                std::cout << "                      (" 
                          << unitCell.lattice_vectors(0,i) << ", " 
                          << unitCell.lattice_vectors(1,i) << ", " 
                          << unitCell.lattice_vectors(2,i) << ") \n";
              }

              continue;
            }
            else if (key == "atom") {
              std::string element;
              lineStream >> element
                         >> unitCell.atomic_positions(0,atom_counter) 
                         >> unitCell.atomic_positions(1,atom_counter) 
                         >> unitCell.atomic_positions(2,atom_counter);
              unitCell.atomic_species.push_back(convertStringToElement(element));
              atom_counter++;
              continue;
            }
            else if (key == "UnitCellDimensions") {
              lineStream >> unitCellDimensions[0] >> unitCellDimensions[1] >> unitCellDimensions[2];
              std::cout << "\n     Unit cell dimensions = " << unitCellDimensions[0] << " x " 
                                                            << unitCellDimensions[1] << " x " 
                                                            << unitCellDimensions[2] << " \n";
              continue;
            }
            else if (key == "InteractionCutoffDistance") {
              lineStream >> interactionCutoffDistance;
              std::cout << "\n     Interaction cutoff distance = " << interactionCutoffDistance << "\n";
              continue;
            }

          }

        }

      }

      inputFile.close();

    }

    // Sanity checks:
    assert(atom_counter == unitCell.number_of_atoms);

    std::cout << "\n     Decorated atoms in a unit cell:\n";
    std::cout << "     Atom          x              y              z           (unit: lattice constant) \n";
    for (unsigned int i=0; i<unitCell.number_of_atoms; i++) {
      printf("     %s    %12.6f   %12.6f   %12.6f\n", convertElementToString(unitCell.atomic_species[i]).c_str(), 
                                                      unitCell.atomic_positions(0,i),
                                                      unitCell.atomic_positions(1,i),
                                                      unitCell.atomic_positions(2,i));

    }

}


void Lattice::constructUnitCellVectors()
{
 
    unitCellVectors.resize(3, numberOfUnitCells);

    unsigned int counter = 0;
    
    for (unsigned int k=0; k<unitCellDimensions[2]; k++) {
      for (unsigned int j=0; j<unitCellDimensions[1]; j++) {
        for (unsigned int i=0; i<unitCellDimensions[0]; i++) {
          unsigned int unitCellIndex = getUnitCellIndex(i, j, k);
          assert (unitCellIndex == counter);
          unitCellVectors(0,unitCellIndex) = int(i);
          unitCellVectors(1,unitCellIndex) = int(j);
          unitCellVectors(2,unitCellIndex) = int(k);
          counter++;
        }
      }
    }

}


void Lattice::initializeAtomicSpecies()
{
  
  unsigned int atomIndex {0};

  globalAtomicSpecies.resize(totalNumberOfAtoms);
    
  for (unsigned int uc=0; uc<numberOfUnitCells; uc++) {
    for (unsigned int atomID=0; atomID<unitCell.number_of_atoms; atomID++) {
      atomIndex = getAtomIndex(uc, atomID);
      globalAtomicSpecies[atomIndex] = unitCell.atomic_species[atomID];
    }
  }
  
  std::cout << "\n   Initialized atomic species in the whole system based on unit cell information. \n";

}


void Lattice::constructGlobalCoordinates()
{
   
  globalAtomicPositions.resize(3, totalNumberOfAtoms);
    
  for (unsigned int uc=0; uc<numberOfUnitCells; uc++) {
    for (unsigned int atomID=0; atomID<unitCell.number_of_atoms; atomID++) {
      unsigned int atomIndex = getAtomIndex(uc, atomID);
      for (unsigned int i=0; i<3; i++)
        globalAtomicPositions(i, atomIndex) = double(unitCellVectors(i,uc)) + unitCell.atomic_positions(i,atomID);

    }
  }
  
  std::cout << "\n   Constructed global coordinates. \n";

}


void Lattice::writeAtomicPositions(const char* filename)
{

  FILE* atomicPositionFile;
  if (filename != NULL)
    atomicPositionFile = fopen(filename, "w");
  else {
    atomicPositionFile = stdout;
    fprintf(atomicPositionFile, "\n   Atomic Positions:\n");
  }

  fprintf(atomicPositionFile, "     Atom      Element       x              y              z\n");
  for (unsigned int atomID = 0; atomID < totalNumberOfAtoms; atomID++)
    fprintf(atomicPositionFile, "   %5d         %2s   %12.6f   %12.6f   %12.6f\n", 
            atomID, convertElementToString(globalAtomicSpecies[atomID]).c_str(), 
            globalAtomicPositions(0, atomID), globalAtomicPositions(1, atomID), globalAtomicPositions(2, atomID));

  if (filename != NULL) fclose(atomicPositionFile);


}


// Note: Assume the reference unit cell is (0,0,0)
void Lattice::constructRelativeCoordinates()
{

    totalNumberOfNeighboringAtoms = unsigned(nearestNeighborUnitCellList.size()) * unitCell.number_of_atoms;
    relativeAtomicPositions.resize(3, totalNumberOfNeighboringAtoms);

    unsigned int counter {0};
    for (unsigned int uc=0; uc<nearestNeighborUnitCellList.size(); uc++) {
      for (unsigned int atomID=0; atomID<unitCell.number_of_atoms; atomID++) {
        unsigned int atomIndex = uc * unitCell.number_of_atoms + atomID;
        //std::cout << "Relative coordinates " << atomIndex << " = ( ";
        for (unsigned int i=0; i<3; i++) {
          assert (atomIndex == counter);
          relativeAtomicPositions(i, atomIndex) = double(relativeUnitCellVectors(i,uc)) + unitCell.atomic_positions(i,atomID);
          //std::cout << relativeAtomicPositions(i, atomIndex) << " ";
        }
        //std::cout << ") \n";
        counter++;
      }
    }

    std::cout << "\n   Constructed relative coordinates. \n";

}


double Lattice::getPairwiseDistance(unsigned int atom1, unsigned int atom2)
{  

    double distance = 0.0;
    for (unsigned int i=0; i<3; i++)
      distance += (globalAtomicPositions(i, atom2) - globalAtomicPositions(i, atom1)) * (globalAtomicPositions(i, atom2) - globalAtomicPositions(i, atom1));

    return sqrt(distance);  

}  
  

double Lattice::getRelativePairwiseDistance(unsigned int atom1, unsigned int atom2)
{  

    double distance = 0.0;
    for (unsigned int i=0; i<3; i++)
      distance += (relativeAtomicPositions(i, atom2) - relativeAtomicPositions(i, atom1)) * (relativeAtomicPositions(i, atom2) - relativeAtomicPositions(i, atom1));

    return sqrt(distance);

  } 


  void Lattice::printAllPairwiseDistances()
  {  
 
    for (unsigned int atom1=0; atom1<totalNumberOfAtoms; atom1++) {
      for (unsigned int atom2=atom1; atom2<totalNumberOfAtoms; atom2++) {
        std::cout << "atom " << atom1 << ", atom " << atom2 << " : " 
                  << globalAtomicPositions(0, atom2) - globalAtomicPositions(0, atom1) << " , " 
                  << globalAtomicPositions(1, atom2) - globalAtomicPositions(1, atom1) << " , " 
                  << globalAtomicPositions(2, atom2) - globalAtomicPositions(2, atom1) << " . " 
                  << getPairwiseDistance(atom1,atom2) << "\n";
      }
    }

  }


  void Lattice::printPairwiseDistancesInUnitCellList(unsigned int atom1_global)
  {

    // Current unit cell is always (0,0,0) relatively.
    unsigned int localUnitCellIndex = getRelativeUnitCellIndex(0, 0, 0);
    unsigned int globalUnitCellIndex = atom1_global / unitCell.number_of_atoms;
    
    unsigned int atom1 = getAtomIndex(localUnitCellIndex, atom1_global % unitCell.number_of_atoms);

    for (unsigned int j=0; j<nearestNeighborUnitCellList[globalUnitCellIndex].size(); j++) {
      for (unsigned int k=0; k<unitCell.number_of_atoms; k++) {
        unsigned int atom2 = getAtomIndex(j, k);
        unsigned int atom2_global = getAtomIndex(nearestNeighborUnitCellList[globalUnitCellIndex][j], k);
        std::cout << "atom " << atom1_global << ", atom " << atom2_global << " : " 
                  << relativeAtomicPositions(0, atom2) - relativeAtomicPositions(0, atom1) << " , " 
                  << relativeAtomicPositions(1, atom2) - relativeAtomicPositions(1, atom1) << " , " 
                  << relativeAtomicPositions(2, atom2) - relativeAtomicPositions(2, atom1) << " . " 
                  << getRelativePairwiseDistance(atom1, atom2) << "\n";
      }
    }

  }


  void Lattice::constructRelativeUnitCellVectors()
  {
    unsigned long temp =  2 * numAdjacentUnitCells + 1;
    unsigned long numberOfNeighoringUnitCells = temp * temp * temp;    // 3 dimensions, including own unit cell.
    relativeUnitCellVectors.resize(3, numberOfNeighoringUnitCells);

    unsigned long counter {0};
    for (int k = -int(numAdjacentUnitCells); k <= int(numAdjacentUnitCells); k++) {
      for (int j = -int(numAdjacentUnitCells); j <= int(numAdjacentUnitCells); j++) {
        for (int i = -int(numAdjacentUnitCells); i <= int(numAdjacentUnitCells); i++) {
          relativeUnitCellVectors(0, counter) = i;
          relativeUnitCellVectors(1, counter) = j;
          relativeUnitCellVectors(2, counter) = k;
          counter++;
        }
      }
    }

    assert(counter == numberOfNeighoringUnitCells);

  }


  std::vector<unsigned int> Lattice::constructNearestNeighborUnitCellList(unsigned int currentUnitCell)
  {

    // A neighbor list where pairwise interactions with current unit cell will be checked
    std::vector<unsigned int> unitCellList;

    unsigned int nx_new, ny_new, nz_new;
    int nx = unitCellVectors(0, currentUnitCell);
    int ny = unitCellVectors(1, currentUnitCell);
    int nz = unitCellVectors(2, currentUnitCell);

    // A lambda to calculate neighboring unit cell components considering p.b.c. shift
    auto getUnitCellComponentPBC = [=](int x_new, unsigned int dimension) { 
      return x_new >= 0 ? (unsigned(x_new) % unitCellDimensions[dimension]) : (unsigned(x_new + int(unitCellDimensions[dimension])));
    };

    for (int k = -int(numAdjacentUnitCells); k <= int(numAdjacentUnitCells); k++) {
      for (int j = -int(numAdjacentUnitCells); j <= int(numAdjacentUnitCells); j++) {
        for (int i = - int(numAdjacentUnitCells); i <= int(numAdjacentUnitCells); i++) {
/*
    // Note: The following 4 lines were used to reduce the number of neighboring unit cells by half. 
    for (int k=-numAdjacentUnitCells; k<=0; k++) {
      for (int j=-numAdjacentUnitCells; j<=numAdjacentUnitCells; j++) {
        for (int i=-numAdjacentUnitCells; i<=numAdjacentUnitCells; i++) {
          if (getRelativeUnitCellIndex(i,j,k) > 13) break;
*/
          nx_new = getUnitCellComponentPBC(nx+i, 0);
          ny_new = getUnitCellComponentPBC(ny+j, 1);
          nz_new = getUnitCellComponentPBC(nz+k, 2);
          //std::cout << i << " " << j << " " << k << " : " << getUnitCellIndex(nx_new, ny_new, nz_new) << "\n";
          unitCellList.push_back(getUnitCellIndex(nx_new, ny_new, nz_new));
        }
      }
    }
    
    return unitCellList;

  }


// Construct the neighbor list for each atom, given the list of neighboring unit cells to check
std::vector<AtomBase> Lattice::constructNeighborListFromNeighboringUnitCells(unsigned int atom1_global)
{

    // A neighbor list where pairwise interactions with current atom will be checked
    std::vector<AtomBase> atomList;

    // Current unit cell is always (0,0,0) relatively.
    unsigned int localUnitCellIndex = getRelativeUnitCellIndex(0, 0, 0);
    unsigned int globalUnitCellIndex =  atom1_global / unitCell.number_of_atoms;
    
    unsigned int atom1 = getAtomIndex(localUnitCellIndex,  atom1_global % unitCell.number_of_atoms);

    for (unsigned int j=0; j<nearestNeighborUnitCellList[globalUnitCellIndex].size(); j++) {
      for (unsigned int k=0; k<unitCell.number_of_atoms; k++) {

        unsigned int atom2 = getAtomIndex(j, k);
        unsigned int atom2_global = getAtomIndex(nearestNeighborUnitCellList[globalUnitCellIndex][j], k);
        if (atom1_global == atom2_global) break;
        double distance = getRelativePairwiseDistance(atom1, atom2);

        if (distance <= interactionCutoffDistance)
          atomList.push_back({atom2_global, distance});
        
      }
    }

    std::sort(atomList.begin(), atomList.end(), 
              [](const auto& a, const auto& b) { return a.atomID < b.atomID; }
    );

    std::cout << "Sorted neighbor list of " << atom1_global << ":\n";
    std::cout << "Atom    distance \n";
    for (auto i : atomList)
      std::cout << i.atomID << " " << i.distance << "\n";

  return atomList;

}


// Construct the neighbor list for each atom in a unit cell
// Unit cell and atomID are relative to the reference unit cell
void Lattice::constructPrimaryNeighborList()
{

  double distance {0.0};
  double dx       {0.0};
  double dy       {0.0};
  double dz       {0.0};

  std::cout << "\n   Construct primary neighbor list:\n";

  primaryNeighborList.resize(unitCell.number_of_atoms);
  unsigned int localUnitCellIndex = getRelativeUnitCellIndex(0, 0, 0);

  for (unsigned int atomID=0; atomID<unitCell.number_of_atoms; atomID++) {
    unsigned int atom1 = getAtomIndex(localUnitCellIndex,  atomID);

    for (unsigned int j=0; j<nearestNeighborUnitCellList[localUnitCellIndex].size(); j++) {
      for (unsigned int k=0; k<unitCell.number_of_atoms; k++) {
        unsigned int atom2 = getAtomIndex(j, k);
        //if (atom1 == atom2) break;                // avoids double counting within the same unit cell
        if (atom1 == atom2) continue;               // avoids putting the reference atom itself into the neighbor list
        dx = relativeAtomicPositions(0, atom2) - relativeAtomicPositions(0, atom1);
        dy = relativeAtomicPositions(1, atom2) - relativeAtomicPositions(1, atom1);
        dz = relativeAtomicPositions(2, atom2) - relativeAtomicPositions(2, atom1);
        distance = getRelativePairwiseDistance(atom1, atom2);

        // Store the distance if it is not yet in neighborDistances
        if (!isFoundInVector(distance, neighborDistances))
          neighborDistances.push_back(distance);

        // Add the atom to neighbor list if within cutoff
        if (distance <= interactionCutoffDistance)
          primaryNeighborList[atomID].push_back({atom2, distance});

      }
    }

    // Print the primary neighbor list for the current atom
    std::cout << "\n     Primary neighbor list of " << atom1 << ":\n";
    std::cout << "          Atom      Distance \n";
    for (auto i : primaryNeighborList[atomID])
      printf("     %8d   %12.6f\n", i.atomID, i.distance);

  }

  // Sort the neighbor distance list and print out
  std::sort(neighborDistances.begin(), neighborDistances.end(), 
            [](const auto& a, const auto& b) { return a < b; }
  );
  std::cout << "\n     Neighboring distances: \n";
  for (unsigned int i=0; i<neighborDistances.size(); i++)
    printf("     %2dth neighbor : %12.6f \n", i, neighborDistances[i]);

  std::cout << "\n   Constructed primary neighbor lists for all atoms in a unit cell. \n";

}


void Lattice::mapPrimaryToAllNeighborLists()
{

  unsigned int thisAtom, atom_tmp, relative_uc, real_uc, atomID_in_uc, atomID;

  for (unsigned int i=0; i<numberOfUnitCells; i++) {
    for (unsigned int j=0; j<unitCell.number_of_atoms; j++) {

      thisAtom = getAtomIndex(i,j);

      for (auto k : primaryNeighborList[j]) {
        atom_tmp = k.atomID;

        relative_uc = atom_tmp / unitCell.number_of_atoms;       // which relative unit cell the neighoring atom in?
        atomID_in_uc = atom_tmp % unitCell.number_of_atoms;      // which atom is the neighoring atom in a unit cell?
        real_uc = nearestNeighborUnitCellList[i][relative_uc];
        atomID = getAtomIndex(real_uc, atomID_in_uc);

        neighborList[thisAtom].push_back({atomID, k.distance});

      }

      std::sort(neighborList[thisAtom].begin(), neighborList[thisAtom].end(), 
                [](const auto& a, const auto& b) { return a.atomID < b.atomID; }
      );

    }
  }

  std::cout << "   Mapped primary neighbor lists to all atoms in system. \n";

}
