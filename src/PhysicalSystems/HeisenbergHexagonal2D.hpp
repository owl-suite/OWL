#ifndef HEISENBERGHEXAGONAL2D_HPP
#define HEISENBERGHEXAGONAL2D_HPP

#include <filesystem>
#include "PhysicalSystemBase.hpp"
#include "Main/Globals.hpp"

class HeisenbergHexagonal2D : public PhysicalSystem {

public :

  HeisenbergHexagonal2D(const char* spinConfigFile = "config_initial.dat", int = 0); 
  ~HeisenbergHexagonal2D();

  //void readCommandLineOptions()                         override;
  void writeConfiguration(int = 0, const char* = NULL)  override;
  void getObservables()                                 override;
  void doMCMove()                                       override;
  void acceptMCMove()                                   override;
  void rejectMCMove()                                   override;

  //void buildMPIConfigurationType()                      override;

private :

  unsigned int Size;

  static const int maxShells = 5;
  int numShells;
  double exchangeParameter[maxShells];
  double uniaxialAnisotropy;
  double externalField[3];

  struct SpinDirection {
    double x;
    double y;
    double z;
  };  

  // Old configuration
  unsigned int CurX, CurY;
  SpinDirection CurType;

  // New configuration
  SpinDirection** spin;          // 2D array because it is a 2D model
  //double spinLength;
  
  bool firstTimeGetMeasures;

  // Private functions
  ObservableType                                                             getExchangeInteractions();
  ObservableType                                                             getExternalFieldEnergy();
  ObservableType                                                             getAnisotropyEnergy();
  std::tuple<ObservableType, ObservableType, ObservableType, ObservableType> getMagnetization();
  ObservableType                                                             getDifferenceInExchangeInteractions();
  ObservableType                                                             getDifferenceInExternalFieldEnergy();
  ObservableType                                                             getDifferenceInAnisotropyEnergy();

  void readHamiltonian(const char* inputFile);
  
  void readSpinConfigFile(const std::filesystem::path& spinConfigFile);
  void initializeSpinConfiguration(int initial);

  ObservableType getShell_1_ExchangeInteractions();
  ObservableType getShell_2_ExchangeInteractions();
  ObservableType getShell_3_ExchangeInteractions();
  ObservableType getShell_4_ExchangeInteractions();
  ObservableType getShell_5_ExchangeInteractions();
  
  ObservableType getShell_1_DifferenceInExchangeInteractions();
  ObservableType getShell_2_DifferenceInExchangeInteractions();
  ObservableType getShell_3_DifferenceInExchangeInteractions();
  ObservableType getShell_4_DifferenceInExchangeInteractions();
  ObservableType getShell_5_DifferenceInExchangeInteractions();
};

#endif
