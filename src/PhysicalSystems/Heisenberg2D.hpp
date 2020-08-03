#ifndef HEISENBERG2D_HPP
#define HEISENBERG2D_HPP

#include <filesystem>
#include "PhysicalSystemBase.hpp"
#include "Main/Globals.hpp"

//class Heisenberg2D : public PhysicalSystem<double> {
class Heisenberg2D : public PhysicalSystem {

public :

  Heisenberg2D(const char* spinConfigFile = "config_initial.dat", int = 0); 
  ~Heisenberg2D();

  //void readCommandLineOptions()                         override;
  void writeConfiguration(int = 0, const char* = NULL)  override;
  void getObservables()                                 override;
  void doMCMove()                                       override;
  void acceptMCMove()                                   override;
  void rejectMCMove()                                   override;

  //void buildMPIConfigurationType()                      override;

private :

  unsigned int Size;
  unsigned int LatticeSize;

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

  void GetMeasuresBruteForce();
  void readSpinConfigFile(const std::filesystem::path& spinConfigFile);
  void initializeSpinConfiguration(int initial);

};

#endif
