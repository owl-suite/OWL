#ifndef CRYSTALBASE_HPP
#define CRYSTALBASE_HPP

#include "Utilities/Matrix.hpp"
#include "Main/Communications.hpp"
#include "Main/Globals.hpp"


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
  
  // Unit cell:
  UnitCell                  unitCell;
  std::vector<unsigned int> unitCellDimensions = {0, 0, 0};  // Number of unit cells in each dimension (nx, ny, nz). Initialize to 0.
  Matrix<int>               unitCellVectors;
  unsigned int              numberOfUnitCells  {0};
  
  // Atoms:
  unsigned int              totalNumberOfAtoms;
  Matrix<double>            globalAtomicPositions;           // Global atomic positions in the crystal (in lattice constant)

  // Neighbor lists (might move to the PhysicalSystem derived class):
  Matrix<int>                              relativeUnitCellVectors;
  std::vector< std::vector<unsigned int> > nearestNeighborUnitCellList;     // Each unit cell has a list of nearest neighbors
  Matrix<double>                           relativeAtomicPositions;         // Relative atomic positions in neighboring unit cells (in lattice constant)
  unsigned int                             totalNumberOfNeighboringAtoms;
  int                                      numAdjacentUnitCells;

  // Constructor 1: initialize unit cell and lattice from input file
  Lattice(const char* inputFile);
  
  // [TODO] Constructor 2 for initializing unit cell in derived classes where crystal structures are known and standard 

  // Destructor
  ~Lattice() { }


  // Member functions:
  
  void readUnitCellInfo(const char* mainInputFile);
  void constructUnitCellVectors();
  void constructGlobalCoordinates();

  inline unsigned int getUnitCellIndex(unsigned int x, unsigned int y, unsigned int z)          // returns a unique ID for a unit cell
  { return z * unitCellDimensions[0] * unitCellDimensions[1] + y * unitCellDimensions[0] + x; }
  
  inline unsigned int getAtomIndex(unsigned int x, unsigned int y, unsigned int z, unsigned int atomIDinUnitCell)
  { return (z * unitCellDimensions[0] * unitCellDimensions[1] + y * unitCellDimensions[0] + x) * unitCell.number_of_atoms + atomIDinUnitCell; }

  inline unsigned int getAtomIndex(unsigned int unitCellIndex, unsigned int atomIDinUnitCell)
  { return unitCellIndex * unitCell.number_of_atoms + atomIDinUnitCell; }

  double getPairwiseDistance(unsigned int atom1, unsigned int atom2);
  void   printAllPairwiseDistances();

  // Neighbor-list related
  std::vector<unsigned int> constructNearestNeighborUnitCellList(unsigned int currentUnitCell);
  void                      constructRelativeUnitCellVectors();
  void                      constructRelativeCoordinates();
  double                    getRelativePairwiseDistance(unsigned int atom1, unsigned int atom2);
  void                      printPairwiseDistancesInUnitCellList(unsigned int atomID);

  inline unsigned int getRelativeUnitCellIndex(unsigned int x, unsigned int y, unsigned int z)
  { 
    unsigned int temp =  2 * numAdjacentUnitCells + 1;
    return (temp * temp * temp / 2) + (z * temp * temp) + (y * temp) + x; 
    //return z * 9 + y * 3 + x + 13;;
  }

};

#endif
