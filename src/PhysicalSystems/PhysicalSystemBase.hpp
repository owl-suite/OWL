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
  virtual void writeConfiguration(int = 0, const char* = NULL) = 0;
  virtual void getObservables() = 0;
  virtual void doMCMove() = 0;
  virtual void acceptMCMove() = 0;
  virtual void rejectMCMove() = 0;        // restore old observables and old configurations to current ones
  
  // Construct data structures for MPI communications. Used in Replica exchanges.
  void buildMPIConfigurationType() {};

  void readCommandLineOptions() {};
  void readHamiltonianTerms(const char* inputFile) {};

  // Parameters common to (needed by) all models:
  unsigned long int numberOfMCSweepsPerStep {1};
  unsigned int numObservables;
  ObservableType* observables;
  ObservableType* oldObservables;

  // For systems where energy is calculated from the difference with the previous configuration, 
  // this flag will cause the energy to be calculated from scratch again.
  // Useful for initialization or after replica exchange.
  bool getObservablesFromScratch {true};   

  // MPI derived type to store configuration for replica exchange
  MPI_Datatype MPI_ConfigurationType;
  void* pointerToConfiguration;

  // MPI Communicator for one energy calculation
  //MPICommunicator PhysicalSystemCommunicator;


protected:

  void initializeObservables(unsigned int n) {
    numObservables = n;
    if (numObservables > 0) {
      observables    = new ObservableType[numObservables];
      oldObservables = new ObservableType[numObservables];
      for (unsigned int i=0; i<numObservables; i++) {
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
    for (unsigned int i=0; i<numObservables; i++) {
      oldObservables[i] = observables[i];
      observables[i] = 0;
    }
  }


  void deleteObservables() {
    delete[] observables;
    delete[] oldObservables;
  }

};

#endif
