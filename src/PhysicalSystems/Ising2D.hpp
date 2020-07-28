#ifndef ISING2D_HPP
#define ISING2D_HPP

#include "PhysicalSystemBase.hpp"
//#include "Main/Globals.hpp"


class Ising2D : public PhysicalSystem {

public :

  Ising2D(const char* = NULL, int = 0); 
  ~Ising2D();

  //void readCommandLineOptions()                         override;
  void writeConfiguration(int = 0, const char* = NULL)  override;
  void getObservables()                                 override;
  void doMCMove()                                       override;
  void acceptMCMove()                                   override;
  void rejectMCMove()                                   override;

  void buildMPIConfigurationType();

private :

  // YingWai: too much trouble to MPI_send an enum... (Oct 10, 2017)
  //enum SpinDirection {DOWN = -1, UP = +1};
  typedef int SpinDirection;

  unsigned int Size;
  unsigned int LatticeSize;

  // Old configuration
  unsigned int CurX, CurY;
  SpinDirection CurType;

  // New configuration
  //SpinDirection** spin;          // 2D array because it is a 2D model
  SpinDirection* spin;             // make it a flat array for MPI to operate on  (Oct 10, 2017)
  

  // YingWai: these do not seem to be needed anymore (Oct 10, 2017)
  //bool firstTimeGetMeasures {true};
  //void GetMeasuresBruteForce();

};

#endif
