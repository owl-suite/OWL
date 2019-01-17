#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <cstdlib>
#include "Main/Globals.hpp"
#include "Main/Communications.hpp"


// To DO: make it a template class to allow for int / double observables  (July 15, 2017)
class PhysicalSystem {

public:
  
  // Constructor:
  PhysicalSystem() {}

  // Destructor:
  virtual ~PhysicalSystem() {}  

  // Functions:
  virtual void readCommandLineOptions() = 0;
  virtual void writeConfiguration(int = 0, const char* = NULL) = 0;
  virtual void getObservables() = 0;
  virtual void doMCMove() = 0;
  virtual void acceptMCMove() = 0;
  virtual void rejectMCMove() = 0;        // restore old observables and old configurations to current ones
  
  virtual void buildMPIConfigurationType() = 0;

  // Parameters common to (needed by) all models:
  int numObservables;
  ObservableType* observables;
  ObservableType* oldObservables;
  //double* observables;
  //double* oldObservables;

  // MPI derived type to store configuration for replica exchange
  MPI_Datatype MPI_ConfigurationType;
  void* pointerToConfiguration;

  // For systems where energy is calculated from the difference with the previous configuration, 
  // this flag will cause the energy to be calculated from scratch again
  // Useful for initialization and after replica exchange
  bool getObservablesFromScratch {true};        

  // MPI Communicator for one energy calculation
  //MPICommunicator PhysicalSystemCommunicator;


protected:

  void initializeObservables(int n) {
    numObservables = n;
    if (numObservables > 0) {
      observables    = new ObservableType[numObservables];
      oldObservables = new ObservableType[numObservables];
      //observables    = new double[numObservables];
      //oldObservables = new double[numObservables];
      for (int i=0; i<numObservables; i++) {
        observables[i]    = 0;
        oldObservables[i] = 0;
      }
    }
    else {
      std::cerr << "Error: Number of observables unphysical." << std::endl;
      exit(5);
    }
  }

  void resetObservables() {
    for (int i=0; i<numObservables; i++) {
      oldObservables[i] = observables[i];
      observables[i] = 0;
    }
  }

/*
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
