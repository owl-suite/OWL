#ifndef GENERAL_3D_CRYSTAL_STRUCTURE_HPP
#define GENERAL_3D_CRYSTAL_STRUCTURE_HPP

#include "PhysicalSystemBase.hpp"
#include "CrystalBase.hpp"
#include "Main/Globals.hpp"


struct SpinDirection {
  double x;
  double y;
  double z;
};


struct NeighboringAtomInfo {
  unsigned int atomID;
  double       distance {0.0};          // Distance from a reference atom
  double       J_ij     {0.0};          // Exchange coupling
  double       D_ij     {0.0};          // Dzyaloshinskii-Moriya (DM) interaction 
};


class CrystalStructure3D : public PhysicalSystem {

public :

  CrystalStructure3D(const char* inputFile, const char* spinConfigFile = NULL, int initial = 0); 
  ~CrystalStructure3D();

  void writeConfiguration(int = 0, const char* = NULL)  override;
  void getObservables()                                 override;
  void doMCMove()                                       override;
  void acceptMCMove()                                   override;
  void rejectMCMove()                                   override;

  //void buildMPIConfigurationType()                      override;

  void readHamiltonianTerms(const char* inputFile);

    Lattice                                         lattice;

private :

  // Overall configuration

  std::vector<SpinDirection>                      spin;
  std::vector< std::vector<NeighboringAtomInfo> > neighborList;            // Each atom has a list of nearest neighbor atoms
  std::vector< std::vector<NeighboringAtomInfo> > primaryNeighborList;     // Neighbor list for each atom in a unit cell 

  // Old configuration
  unsigned int  currentPosition;
  SpinDirection currentSpin;
  
  bool firstTimeGetMeasures;


  // Private functions
  // Initialization:
  void readSpinConfigFile(const char* spinConfigFile);
  void initializeSpinConfiguration(int initial);
  void assignRandomSpinConfiguration(unsigned int currentAtom);
  std::vector<NeighboringAtomInfo> constructNeighborListFromNeighboringUnitCells(unsigned int currentAtom);

  void constructPrimaryNeighborList();
  void mapPrimaryToAllNeighborLists();

  double assignExchangeCouplings(double dx, double dy, double dz, double dr);
  double assignDMInteractions(double dx, double dy, double dz, double dr);

  // Hamiltonian measurements:
  void getObservablesFromScratch();
  ObservableType nearestNeighborInterations();
  
};

#endif
