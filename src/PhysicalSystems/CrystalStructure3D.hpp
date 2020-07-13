#ifndef GENERAL_3D_CRYSTAL_STRUCTURE_HPP
#define GENERAL_3D_CRYSTAL_STRUCTURE_HPP

#include <filesystem>
#include <tuple>
#include "PhysicalSystemBase.hpp"
#include "CrystalBase.hpp"
#include "Main/Globals.hpp"
#include "Utilities/CompareNumbers.hpp"


struct SpinDirection {
  double x;
  double y;
  double z;
};


struct NeighboringAtomInfo {
  unsigned int atomID;
  double       distance {0.0};          // Distance from a reference atom
  double       J_ij     {0.0};          // Exchange coupling
  double       D_ij     {0.0};          // Dzyaloshinskii-Moriya (DM) interaction in z-direction
};


class CrystalStructure3D : public PhysicalSystem {

public :

  CrystalStructure3D(const char* inputFile, const std::filesystem::path& spinConfigFile = std::filesystem::path(), int initial = 0); 
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
  std::vector< std::vector<NeighboringAtomInfo> > neighborList;                       // Each atom has a list of neighboring atoms
  std::vector< std::vector<NeighboringAtomInfo> > primaryNeighborList;                // Neighbor list for each atom in a unit cell 

  std::vector<double>                             neighborDistances;                  // Stores the distances between neighbors
  double                                          interactionCutoffDistance {1.0};    // Default to one lattice constant

  // Old configuration
  unsigned int  currentPosition;
  SpinDirection oldSpin;
  
  bool firstTimeGetMeasures;


  // Private functions

  // Initialization:
  void   readSpinConfigFile(const std::filesystem::path& spinConfigFile);
  void   readInteractionCutoffDistance(const char* inputFile);
  void   initializeSpinConfiguration(int initial);
  void   assignRandomSpinDirection(unsigned int currentAtom);
  double assignExchangeCouplings(double dx, double dy, double dz, double dr);
  //double assignExchangeCouplings_testing(double dx, double dy, double dz, double dr);       // temp. testing code
  double assignDzyaloshinskiiMoriyaInteractions(double dx, double dy, double dz, double dr);
  
  // Neighbor lists:
  std::vector<NeighboringAtomInfo> constructNeighborListFromNeighboringUnitCells(unsigned int currentAtom);
  void                             constructPrimaryNeighborList();
  void                             mapPrimaryToAllNeighborLists();

  // Hamiltonian measurements:
  void                                                       getObservablesFromScratch();
  ObservableType                                             getExchangeInterations();
  ObservableType                                             getDzyaloshinskiiMoriyaInterations();
  std::tuple<ObservableType, ObservableType, ObservableType> getMagnetization();

  ObservableType                                             getDifferenceInExchangeInterations();
  ObservableType                                             getDifferenceInDzyaloshinskiiMoriyaInterations();
};

#endif
