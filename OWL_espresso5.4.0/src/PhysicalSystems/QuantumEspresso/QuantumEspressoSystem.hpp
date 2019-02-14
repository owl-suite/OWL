#ifndef QUANTUM_ESPRESSO_SYSTEM_HPP
#define QUANTUM_ESPRESSO_SYSTEM_HPP

#include <mpi.h>
#include <vector>
#include "PhysicalSystems/PhysicalSystemBase.hpp"
#include "Utilities/Matrix.hpp"
#include "Main/Globals.hpp"


// YingWai's Note:  (Dec 26, 17)
// If this is changed, the MPI derived type (MPIConfigurationType) built by buildMPIConfigurationType() needs to be modified too.
struct QEConfiguration
{
  Matrix<double> atomic_positions;            // Atomic positions (in Angstrom)
  Matrix<double> lattice_vectors;             // Unit cell vectors (in Angstrom)
  std::vector<int> atomic_species;            // Atomic species (in atomic number)
};


class QuantumEspressoSystem : public PhysicalSystem {

public:

  // Constructor
  QuantumEspressoSystem(MPICommunicator PhysicalSystemComm);
  // Destructor
  ~QuantumEspressoSystem();

  void readCommandLineOptions();
  void writeConfiguration(int = 0, const char* = NULL);
  void getObservables();
  void doMCMove();
  void acceptMCMove();
  void rejectMCMove();

  void buildMPIConfigurationType();

private:

  // Parameters from command line options for controlling parallelization in Quantum Espresso
  int   nimage            {1};                 // number of images
  int   npool             {1};                 // (?)
  int   ndiag             {1};                 // (?)
  int   ntg               {1};                 // ntaskg (?)
  int   nband             {1};                 // number of bands
  char  QEInputFile[81]   = { ' ', '\0' };     // QE input file

  // System information
  int natom;                                   // Total number of atoms in system
                                               // (total energy measured in Ry)
  QEConfiguration trialConfig;                 // Trial configuration
  QEConfiguration oldConfig;                   // Old configuration for restoration after a rejected MC move

  // I/O
  void writeSystemFile(const char* = NULL);
  void writeQErestartFile(const char* = NULL);


  // Temp. parameters for MPI
  int MPI_exit_status;

};

#endif

