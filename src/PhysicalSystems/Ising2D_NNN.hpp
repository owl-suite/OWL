#ifndef ISING2D_NNN_HPP
#define ISING2D_NNN_HPP

#include <filesystem>
#include "PhysicalSystemBase.hpp"

class Ising2D_NNN : public PhysicalSystem {

public :

  Ising2D_NNN(const char* spinConfigFile = "config_initial.dat"); 
  ~Ising2D_NNN();

  //void readCommandLineOptions()                         override;
  void writeConfiguration(int = 0, const char* = NULL)  override;
  void getObservables()                                 override;
  void doMCMove()                                       override;
  void acceptMCMove()                                   override;
  void rejectMCMove()                                   override;

  void buildMPIConfigurationType();

private :

  unsigned int                numExchangeInteractions {2};         // J1 and J2
  std::vector<ObservableType> exchangeInteraction;

  typedef int SpinDirection;

  unsigned int Size;

  // Old configuration
  unsigned int CurX, CurY;
  SpinDirection oldSpin;

  // New configuration
  SpinDirection* spin;             // make it a flat array for MPI to operate on

  // Initialization:
  void   readSpinConfigFile(const std::filesystem::path& spinConfigFile);
  void   readHamiltonian(const std::filesystem::path& mainInputFile);
  void   initializeSpinConfiguration(unsigned int initial = 0);

};

#endif
