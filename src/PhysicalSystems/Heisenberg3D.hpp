#ifndef HEISENBERG3D_HPP
#define HEISENBERG3D_HPP

#include <filesystem>
#include "PhysicalSystemBase.hpp"
#include "Main/Globals.hpp"

class Heisenberg3D : public PhysicalSystem {

public :

  Heisenberg3D(const char* spinConfigFile = "config_initial.dat", int initial = 0); 
  ~Heisenberg3D();

  //void readCommandLineOptions()                         override;
  void writeConfiguration(int = 0, const char* = NULL)  override;
  void getObservables()                                 override;
  void doMCMove()                                       override;
  void acceptMCMove()                                   override;
  void rejectMCMove()                                   override;

  //void buildMPIConfigurationType()                      override;

private :

  unsigned int Size;

  struct SpinDirection {
    double x;
    double y;
    double z;
  };  
 
  // Interaction constants (default to AFM Heisenberg model nearest-neighbor interactions)
  InteractionConstant J1 {-1.0};
  InteractionConstant J2 {-1.0};
  InteractionConstant J3 {0.0};
  InteractionConstant J4 {0.0};
  InteractionConstant Dz {0.0};
  
  // Old configuration
  unsigned int CurX, CurY, CurZ;
  SpinDirection oldSpin;

  // New configuration
  SpinDirection*** spin;          // 3D array because it is a 3D model
  //double spinLength;
  
  bool firstTimeGetMeasures;
  
  // Private functions
  ObservableType                                                             getExchangeInteractions();
  ObservableType                                                             getExternalFieldEnergy();
  std::tuple<ObservableType, ObservableType, ObservableType, ObservableType> getMagnetization();
  ObservableType                                                             getDifferenceInExchangeInteractions();
  ObservableType                                                             getDifferenceInExternalFieldEnergy();

  void readInteractionsFromFile(const std::filesystem::path& interactionFile);
  void readSpinConfigFile(const std::filesystem::path& spinConfigFile);
  void initializeSpinConfiguration(int initial);

};

#endif
