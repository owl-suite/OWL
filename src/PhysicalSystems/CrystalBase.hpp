#ifndef CRYSTALBASE_HPP
#define CRYSTALBASE_HPP

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>
#include "Utilities/Matrix.hpp"
#include "Main/Communications.hpp"
#include "Main/Globals.hpp"


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



// Note: This is the same as the QEConfiguration struct in QuantumEspressoSystem.hpp
struct UnitCell
{
  Matrix<double>           lattice_vectors;              // Unit cell vectors (in Angstrom)
  unsigned int             number_of_atoms;              // Number of atoms in a unit cell
  std::vector<std::string> atomic_species;               // Atomic species
  Matrix<double>           atomic_positions;             // Atomic positions (in lattice constant)
};


class Lattice {

public :
  
  // Member objects:
  
  UnitCell                  unitCell;
  std::vector<unsigned int> unitCellDimensions = {0, 0, 0};   // number of unit cells in each dimension (nx, ny, nz). Initialize to 0.
  Matrix<unsigned int>      unitCellVectors;
  unsigned int              numberOfUnitCells  {0};
  unsigned int              unitCellIndex;
  
  Matrix<double>            globalAtomicPositions;        // Global atomic positions in the crystal (in lattice constant)


  // Constructor 1: initialize unit cell and lattice from input file
  Lattice(const char* inputFile)
  {
    
    // Initialize unit cell
    unitCell.lattice_vectors.resize(3, 3);
    readUnitCellInfo(inputFile);

    numberOfUnitCells = unitCellDimensions[0] * unitCellDimensions[1] * unitCellDimensions[2];
    constructUnitCellVectors();

    // Initialize crystal
    constructGlobalCoordinates();

    
  }


  // TODO: Constructor for initializing unit cell in derived classes where crystal structures are known and standard 

  // Destructor
  ~Lattice() { }


  // Member functions:
  
  // Note: unit cell layout in lattice_vectors(i,j):
  //       j ->
  //     i    ax  bx  cx
  //     |    ay  by  cy
  //     v    az  bz  cz
  void readUnitCellInfo(const char* mainInputFile)
  {

    //if (GlobalComm.thisMPIrank == 0)
      if (std::filesystem::exists(mainInputFile))
        std::cout << "Crystal class reading input file: " << mainInputFile << "\n";

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

            if (key == "numberOfAtomsPerUnitCell") {
              lineStream >> unitCell.number_of_atoms;
              std::cout << "Crystal: number of atoms per unit cell = " << unitCell.number_of_atoms << "\n";

              // Allocate memory for atomic_positions
              unitCell.atomic_positions.resize(3, unitCell.number_of_atoms);
              continue;
            }
            else if (key == "latticeVectors") {
              lineStream >> unitCell.lattice_vectors(0,0) >> unitCell.lattice_vectors(1,0) >> unitCell.lattice_vectors(2,0);
              std::cout << "Crystal: Lattice vectors (" 
                        << unitCell.lattice_vectors(0,0) << ", " 
                        << unitCell.lattice_vectors(1,0) << ", " 
                        << unitCell.lattice_vectors(2,0) << ") \n";

              // Read the other two lattice vectors
              for (unsigned int i=1; i<3; i++) {
                lineStream.clear();
                std::getline(inputFile, line);               
                if (!line.empty()) lineStream.str(line);
                lineStream >> unitCell.lattice_vectors(0,i) >> unitCell.lattice_vectors(1,i) >> unitCell.lattice_vectors(2,i);
                std::cout << "Crystal: Lattice vectors (" 
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
              //std::cout << "element = " << element << "\n";
              unitCell.atomic_species.push_back(element);
              atom_counter++;
              //std::cout << "atom_counter = " << atom_counter << "\n";
              continue;
            }
            else if (key == "unitCellDimensions") {
              lineStream >> unitCellDimensions[0] >> unitCellDimensions[1] >> unitCellDimensions[2];
              std::cout << "Crystal: unit cell dimensions = " << unitCellDimensions[0] << " x " 
                                                              << unitCellDimensions[1] << " x " 
                                                              << unitCellDimensions[2] << " \n";
              continue;
            }

          }

        }

      }

      inputFile.close();

    }

    // Sanity checks:
    assert(atom_counter == unitCell.number_of_atoms);

    for (unsigned int i=0; i<unitCell.number_of_atoms; i++) {
      std::cout << unitCell.atomic_species[i] << " " 
                << unitCell.atomic_positions(0,i) << " " 
                << unitCell.atomic_positions(1,i) << " " 
                << unitCell.atomic_positions(2,i) << "\n";
    }

  }


  // returns a unique ID for a unit cell
  inline unsigned int getUnitCellIndex(unsigned int x, unsigned int y, unsigned int z)
  {
    return int(z * unitCellDimensions[0] * unitCellDimensions[1] + y * unitCellDimensions[0] + x);
  }


  void constructUnitCellVectors()
  {
 
    unitCellVectors.resize(3, numberOfUnitCells);

    unsigned int counter = 0;
    
    for (unsigned int k=0; k<unitCellDimensions[2]; k++) {
      for (unsigned int j=0; j<unitCellDimensions[1]; j++) {
        for (unsigned int i=0; i<unitCellDimensions[0]; i++) {
          unitCellIndex = getUnitCellIndex(i, j, k);
          assert (unitCellIndex == counter);
          unitCellVectors(0,unitCellIndex) = i;
          unitCellVectors(1,unitCellIndex) = j;
          unitCellVectors(2,unitCellIndex) = k;
          counter++;
        }
      }
    }

  }


  void constructGlobalCoordinates()
  {
   
    globalAtomicPositions.resize(3, numberOfUnitCells * unitCell.number_of_atoms);
    
    for (unsigned int k=0; k<numberOfUnitCells; k++) {
      for (unsigned int j=0; j<unitCell.number_of_atoms; j++) {
        unsigned int index = unitCell.number_of_atoms * k + j;
        //std::cout << "Global coordinates " << k << " " << j << " = ( ";
        std::cout << "Global coordinates " << index << " = ( ";
        for (unsigned int i=0; i<3; i++) {
          globalAtomicPositions(i,index) = unitCellVectors(i,k) + unitCell.atomic_positions(i,j);
          std::cout << globalAtomicPositions(i,index) << " ";
        }
        std::cout << ") \n";
      }
    }
  
  }

};




#endif