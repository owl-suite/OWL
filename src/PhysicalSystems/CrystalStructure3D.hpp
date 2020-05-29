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

private :

  Lattice lattice;

  // Old configuration
  unsigned int currentPosition;
  SpinDirection currentSpin;

  // New configuration
  SpinDirection* spin;
  //double spinLength;
  
  bool firstTimeGetMeasures;

  // Private functions
  void readSpinConfigFile(const char* spinConfigFile);
  void initializeSpinConfiguration(int initial);
  void assignRandomSpinConfiguration(unsigned int i);
  void GetMeasuresBruteForce();

  // Hamiltonian Terms
  ObservableType NearestNeighborInterations();
  
};

#endif
