#ifndef ISING2D_HPP
#define ISING2D_HPP

#include <filesystem>
#include "PhysicalSystemBase.hpp"

class Ising2D : public PhysicalSystem {

public :

  Ising2D(const char* spinConfigFile = "config_initial.dat", int = 0); 
  ~Ising2D();

  //void readCommandLineOptions()                         override;
  void writeConfiguration(int = 0, const char* = NULL)  override;
  void getObservables()                                 override;
  void doMCMove()                                       override;
  void acceptMCMove()                                   override;
  void rejectMCMove()                                   override;

  void buildMPIConfigurationType();

private :

  typedef int SpinDirection;

  unsigned int Size;

  // Old configuration
  unsigned int CurX, CurY;
  SpinDirection CurType;

  // New configuration
  SpinDirection* spin;             // make it a flat array for MPI to operate on

  // Initialization:
  void   readSpinConfigFile(const std::filesystem::path& spinConfigFile);
  void   initializeSpinConfiguration(int initial = 0);

};

#endif
