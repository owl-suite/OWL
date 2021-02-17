#ifndef CRYSTALBASE_HPP
#define CRYSTALBASE_HPP

#include "Elements.hpp"
#include "Utilities/Matrix.hpp"
#include "Main/Communications.hpp"
#include "Main/Globals.hpp"


// Note: This is the same as the QEConfiguration struct in QuantumEspressoSystem.hpp
struct UnitCell
{
  Matrix<double>           lattice_vectors;              // Unit cell vectors (in Angstrom)
  unsigned int             number_of_atoms;              // Number of atoms in a unit cell
  std::vector<Element>     atomic_species;               // Atomic species
  Matrix<double>           atomic_positions;             // Atomic positions (in lattice constant)
};


struct AtomBase {
  unsigned int atomID;
  double       distance {0.0};                           // Distance from a reference atom
};


class Lattice {

public :
  
  // Unit cell:
  UnitCell                  unitCell;
  std::vector<unsigned int> unitCellDimensions = {0, 0, 0};  // Number of unit cells in each dimension (nx, ny, nz). Initialize to 0.
  unsigned int              numberOfUnitCells  {0};
  Matrix<int>               unitCellVectors;
  Matrix<int>               relativeUnitCellVectors;

  // Atoms:
  unsigned int              totalNumberOfAtoms;
  Matrix<double>            globalAtomicPositions;           // Global atomic positions in the crystal (in lattice constant)
  std::vector<Element>      globalAtomicSpecies;

  // Neighbor lists:
  std::vector< std::vector<unsigned int> > nearestNeighborUnitCellList;     // Each unit cell has a list of nearest neighbors
  Matrix<double>                           relativeAtomicPositions;         // Relative atomic positions in neighboring unit cells (in lattice constant)
  unsigned int                             totalNumberOfNeighboringAtoms;
  unsigned int                             numAdjacentUnitCells;

  std::vector< std::vector<AtomBase> >     primaryNeighborList;             // Neighbor list for each atom in a unit cell 
  std::vector< std::vector<AtomBase> >     neighborList;                    // Each atom has a list of neighboring atoms
  std::vector<double>                      neighborDistances;               // Stores the distances between neighbors
  std::vector<unsigned int>                coordinationNumbers;

  // [TODO]: 
  // 1. interactionCutoffDistance should be incorporated into the constructor of the Hamiltonian class later when it is implemented,
  //    together with the reading of Hamiltonian terms. (July 7, 20)
  // 2. set cutoff distance to nearest-neighbor only by default
  double interactionCutoffDistance {1.0};                                   // Default to one lattice constant

  // Constructor 1: initialize unit cell and lattice from input file
  Lattice(const char* inputFile);
  
  // [TODO] Constructor 2 for initializing unit cell in derived classes where crystal structures are known and standard 

  // Destructor
  ~Lattice() { }


  // Member functions:
  
  void readUnitCellInfo(const char* mainInputFile);
  void constructUnitCellVectors();
  void initializeAtomicSpecies();
  void constructGlobalCoordinates();
  void writeAtomicPositions(const char* filename = NULL);

  // returns a unique ID of a unit cell (a.k.a. unit cell index) from the unit cell coordinates
  inline unsigned int getUnitCellIndex(unsigned int x, unsigned int y, unsigned int z)
  { return z * unitCellDimensions[0] * unitCellDimensions[1] + y * unitCellDimensions[0] + x; }
  
  // Get atomID in the lattice from the unit cell coordinates and the atom's ID in a unit cell
  inline unsigned int getAtomIndex(unsigned int x, unsigned int y, unsigned int z, unsigned int atomIDinUnitCell)
  { return (z * unitCellDimensions[0] * unitCellDimensions[1] + y * unitCellDimensions[0] + x) * unitCell.number_of_atoms + atomIDinUnitCell; }

  // Get atomID in the lattice from the unit cell index and the atom's ID in a unit cell
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

  std::vector<AtomBase> constructNeighborListFromNeighboringUnitCells(unsigned int currentAtom);
  void                  constructPrimaryNeighborList();
  void                  mapPrimaryToAllNeighborLists();
  void                  getCoordinationNumbers();

  inline unsigned int getRelativeUnitCellIndex(unsigned int x, unsigned int y, unsigned int z)
  { 
    unsigned int temp =  2 * numAdjacentUnitCells + 1;
    return (temp * temp * temp / 2) + (z * temp * temp) + (y * temp) + x; 
    //return z * 9 + y * 3 + x + 13;;
  }

};

#endif
