#ifndef ISING2D_HPP
#define ISING2D_HPP

#include "PhysicalSystemBase.hpp"
#include "Globals.hpp"

//class Ising2D : public PhysicalSystem<int> {
class Ising2D : public PhysicalSystem {

public :

  Ising2D(const char* = NULL, int = 0); 
  ~Ising2D();

  void readCommandLineOptions();
  void writeConfiguration(int = 0, const char* = NULL);
  void getObservables();
  void doMCMove();
  //void undoMCMove();
  void acceptMCMove();
  void rejectMCMove();


private :

  enum SpinDirection {DOWN = -1, UP = +1};

  int Size;
  long LatticeSize;


  // Old configuration
  int CurX, CurY;
  SpinDirection CurType;

  // New configuration
  SpinDirection** spin;          // 2D array because it is a 2D model
  
  bool firstTimeGetMeasures;

  void GetMeasuresBruteForce();

};

#endif
