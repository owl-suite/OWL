#ifndef GENERAL_3D_CRYSTAL_STRUCTURE_HPP
#define GENERAL_3D_CRYSTAL_STRUCTURE_HPP

#include <filesystem>
#include <tuple>
#include <vector>
#include "CrystalBase.hpp"
#include "PhysicalSystemBase.hpp"
#include "Main/Globals.hpp"

struct SpinDirection {
  double x {0.0};
  double y {0.0};
  double z {0.0};
};


struct NeighboringAtom {
  unsigned int atomID;
  double       distance {0.0};                           // Distance from a reference atom
  double       J_ij     {0.0};                           // Exchange coupling
  double       D_ij     {0.0};                           // Dzyaloshinskii-Moriya (DM) interaction in z-direction
};


class CrystalStructure3D : public PhysicalSystem {

public :

  CrystalStructure3D(const char* inputFile, int initial = 0); 
  ~CrystalStructure3D();

  void writeConfiguration(int = 0, const char* = NULL)  override;
  void getObservables()                                 override;
  void doMCMove()                                       override;
  void acceptMCMove()                                   override;
  void rejectMCMove()                                   override;

  void getAdditionalObservables()                       override;

  //void buildMPIConfigurationType()                      override;

  //void readHamiltonianTerms(const char* inputFile);

  Lattice                                         lattice;

private :

  // Model specific information add onto neighborList
  std::vector< std::vector<NeighboringAtom> > primaryNeighborList;
  std::vector< std::vector<NeighboringAtom> > neighborList;

  // Overall configuration
  std::vector<SpinDirection>                  spin;
  std::vector<ObservableType>                 localWindingNumber;

  // Old configuration
  unsigned int  currentPosition;
  SpinDirection oldSpin;
  
  bool firstTimeGetMeasures;


  // Private functions

  // Initialization:
  void   readSpinConfigFile(const std::filesystem::path& spinConfigFile);
  void   initializeSpinConfiguration(int initial = 0);
  void   assignRandomSpinDirection(unsigned int currentAtom);
  double assignExchangeCouplings(double dx, double dy, double dz, double dr);
  //double assignExchangeCouplings_testing(double dx, double dy, double dz, double dr);       // temp. testing code
  double assignDzyaloshinskiiMoriyaInteractions(double dz, double dr);
  
  void   addInteractionsToPrimaryNeighborList();
  void   mapPrimaryToAllNeighborLists();           // TODO: this should override the same function in CrystalBase class

  // Hamiltonian measurements:
  void                                                                       getObservablesFromScratch();
  ObservableType                                                             getExchangeInteractions();
  ObservableType                                                             getDzyaloshinskiiMoriyaInteractions();
  std::tuple<ObservableType, ObservableType, ObservableType, ObservableType> getMagnetization();
  ObservableType                                                             getTotalWindingNumber();

  ObservableType                                                             getDifferenceInExchangeInteractions();
  ObservableType                                                             getDifferenceInDzyaloshinskiiMoriyaInteractions();
  ObservableType                                                             getDifferenceInWindingNumber();
  ObservableType                                                             calculateLocalWindingNumber(unsigned int atomID);

};

#endif
