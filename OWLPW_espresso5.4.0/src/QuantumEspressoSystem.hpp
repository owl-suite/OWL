#ifndef QUANTUM_ESPRESSO_SYSTEM_HPP
#define QUANTUM_ESPRESSO_SYSTEM_HPP

#include "PhysicalSystemBase.hpp"
#include "Matrix.hpp"
#include "Globals.hpp"

class QuantumEspressoSystem : public PhysicalSystem {

public:

  QuantumEspressoSystem(SimulationInfo& sim_info);
  ~QuantumEspressoSystem();

  void readCommandLineOptions(SimulationInfo& sim_info);
  void writeConfiguration(int = 0, const char* = NULL);
  void getObservables();
  void doMCMove();
  //void undoMCMove();
  void acceptMCMove();
  void rejectMCMove();


private:

  // Parameters from command line options for controlling parallelization in Quantum Espresso
  int   nimage            {1};                 // number of images
  int   npool             {1};                 // (?)
  int   ndiag             {1};                 // (?)
  int   ntg               {1};                 // ntaskg (?)
  int   nband             {1};                 // number of bands
  char  QEInputFile[81]   = { ' ', '\0' };     // QE input file

  // System information
  int natom;                       // Total number of atoms in system

  // The following two might not be needed, since they should be the same as observables[0] and trialObservables[0]
  double trialEnergy;              // Trial total energy of system (in Ry)
  double oldEnergy;                // Old total energy of the system (in Ry)

  Matrix<double> trialPos;         // Trial atom position    (in angstrom)
  Matrix<double> oldPos;           // Old atom position for restoration after a rejected MC move      (in angstrom)

  Matrix<double> trialLatticeVec;  // Trial cell vector   (in angstrom)
  Matrix<double> oldLatticeVec;    // Old cell vector for restoration after a rejected MC move (i    n angstrom)

  void writeSystemFile(const char* = NULL);
  void writeQErestartFile(const char* = NULL);


  // Temp. parameters for MPI
  int MPI_exit_status;

};

#endif

