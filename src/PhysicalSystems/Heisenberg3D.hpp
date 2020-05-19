#ifndef HEISENBERG3D_HPP
#define HEISENBERG3D_HPP

#include "PhysicalSystemBase.hpp"
#include "Main/Globals.hpp"

class Heisenberg3D : public PhysicalSystem {

public :

  Heisenberg3D(const char* = NULL, int = 0); 
  ~Heisenberg3D();

  void readCommandLineOptions();
  void writeConfiguration(int = 0, const char* = NULL);
  void getObservables();
  void doMCMove();
  void acceptMCMove();
  void rejectMCMove();

  void buildMPIConfigurationType();

private :

  int Size;
  long LatticeSize;

  struct SpinDirection {
    double x;
    double y;
    double z;
  };  

  // Old configuration
  int CurX, CurY, CurZ;
  SpinDirection CurType;

  // New configuration
  SpinDirection*** spin;          // 3D array because it is a 3D model
  //double spinLength;
  
  bool firstTimeGetMeasures;

  void GetMeasuresBruteForce();

};

#endif
