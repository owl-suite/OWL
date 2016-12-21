#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <cstdlib>

class PhysicalSystem {

public:
  
  // Constructor:
  PhysicalSystem() {}

  // Destructor:
  virtual ~PhysicalSystem() {}  

  // Functions:
  virtual void writeConfiguration(int = 0, const char* = NULL) = 0;
  virtual void getObservables() = 0;
  virtual void doMCMove() = 0;
  //virtual void undoMCMove() = 0;
  virtual void acceptMCMove() = 0;
  virtual void rejectMCMove() = 0;

  // Parameters common to (needed by) all models:
  int numObservables;
  double* observables;
  double* oldObservables;

protected:

  void initializeObservables(int n) {
    numObservables = n;
    if (numObservables > 0) {
      observables    = new double[numObservables];
      oldObservables = new double[numObservables];
      for (int i=0; i<numObservables; i++) {
        observables[i]    = 0.0;
        oldObservables[i] = 0.0;
      }
    }
    else {
      std::cerr << "Error: Number of observables unphysical." << std::endl;
      exit(5);
    }
  }

/*
  void resetObservables() {
    for (int i=0; i<numObservables; i++) {
      oldObservables[i] = observables[i];
      observables[i] = 0.0;
    }
  }

  void restoreObservables() {
    double tmp;
    for (int i=0; i<numObservables; i++) {
      tmp = observables[i];
      observables[i] = oldObservables[i];
      oldObservables[i] = tmp;
    }
  }
*/

  void deleteObservables() {
    delete[] observables;
    delete[] oldObservables;
  }

};

#endif
