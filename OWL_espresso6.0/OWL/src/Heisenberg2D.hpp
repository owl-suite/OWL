#ifndef HEISENBERG2D_HPP
#define HEISENBERG2D_HPP

#include "PhysicalSystemBase.hpp"
#include "Globals.hpp"

class Heisenberg2D : public PhysicalSystem {

public :

  Heisenberg2D(SimulationInfo& sim_info, const char* = NULL, int = 0); 
  ~Heisenberg2D();

  void readCommandLineOptions(SimulationInfo& sim_info);
  void writeConfiguration(int = 0, const char* = NULL);
  void getObservables();
  void doMCMove();
  //void undoMCMove();
  void acceptMCMove();
  void rejectMCMove();


private :

  int Size;
  long LatticeSize;

  struct SpinDirection {
    double x;
    double y;
    double z;
  };  

  // Old configuration
  int CurX, CurY;
  SpinDirection CurType;

  // New configuration
  SpinDirection** spin;          // 2D array because it is a 2D model
  double spinLength;
  
  bool firstTimeGetMeasures;

  void GetMeasuresBruteForce();

};

#endif
