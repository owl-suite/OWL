#ifndef MODEL_HPP
#define MODEL_HPP

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include "Main/Globals.hpp"
#include "Main/Communications.hpp"


// To DO: make it a template class to allow for int / double observables  (July 15, 2017)
class PhysicalSystem {

public:
  
  // Constructor:
  PhysicalSystem() {

    if (!std::filesystem::exists("configurations")) 
      std::filesystem::create_directory("configurations");

  }

  // Destructor:
  virtual ~PhysicalSystem() {}  

  // Functions:
  virtual void writeConfiguration(int = 0, const char* = NULL) = 0;
  virtual void getObservables() = 0;
  virtual void doMCMove() = 0;
  virtual void acceptMCMove() = 0;
  virtual void rejectMCMove() = 0;        // restore old observables and old configurations to current ones
  
  virtual void getAdditionalObservables() {};
  virtual void calculateThermodynamics(std::vector<ObservableType>, std::vector<ObservableType>, double);

  // Construct data structures for MPI communications. Used in Replica exchanges.
  void buildMPIConfigurationType() {};

  //void readHamiltonianTerms(const char* inputFile) {};


  // Parameters common to (needed by) all models:
  unsigned int                systemSize;
  unsigned int                numObservables;
  std::vector<ObservableType> observables;
  std::vector<ObservableType> oldObservables;
  std::vector<std::string>    observableName;

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

  void setSystemSize(unsigned int n) {
    systemSize = n;
  }

  void setSystemSize(unsigned int n, unsigned int d) {
    systemSize = n;
    for (unsigned int i=0; i<d-1; i++)
      systemSize *= n;
  }

  void initializeObservables(unsigned int n) {
    numObservables = n;
    if (numObservables > 0) {
      observables.resize(numObservables, 0.0);
      oldObservables.resize(numObservables, 0.0);
    }
    else {
      std::cerr << "Error: Number of observables unphysical." << "\n";
      exit(5);
    }
  }


  void resetObservables() {
    for (unsigned int i=0; i<numObservables; i++) {
      oldObservables[i] = observables[i];
      observables[i] = 0;
    }
  }

};

#endif
