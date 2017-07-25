#ifndef COMMUNICATIONS_HPP
#define COMMUNICATIONS_HPP

#include "mpi.h"
#include "Globals.hpp"


class MPICommunicator {

public :

  MPI_Comm communicator;
  int      thisMPIrank;
  int      totalMPIranks;

  void initialize(MPI_Comm incomingComm);    // this initializes communicator
  void finalize();                           // this finalizes communicator

};

// Global variables for MPI communicators info
extern MPICommunicator GlobalComm;


// Functions
void initializeMPICommunication(SimulationInfo, MPICommunicator&, MPICommunicator&);
void finalizeMPICommunication(SimulationInfo);

// The following should be replaced once the MPI classes are ready? (Jun 14, 17)
void initializeQEMPICommunication();
void finalizeQEMPICommunication();



#endif
