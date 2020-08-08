//  This header file defines a struct "SimulationInfo" that stores all the global variables / data structures
//  for sharing among different part of the code

#ifndef GLOBALS_HPP
#define GLOBALS_HPP

/// Global data types
//enum type_name {integer = 1, single_precision = 2, double_precision = 3};
//typedef int ObservableType;

typedef double       ObservableType;
typedef double       InteractionConstant;
typedef unsigned int indexType;
typedef double       InteractionConstant;


struct SimulationInfo
{

  // Basic simulation info
  int  restartFlag           {0};
  int  algorithm             {2};
  int  system                {2};
  int  rngSeed               {-1};
  char* MCInputFile;
  char* HistogramCheckpointFile;

  // MPI info
  int  numWalkers            {1};
  int  numMPIranksPerWalker  {1};
  int  myWalkerID            {-1};
  int  numberOfWindows       {-1};                   // only valid for REWL, will remain -1 for other MC algorithms
  int  numberOfWalkersPerWindow {-1};                // ditto

  // Command line options for QE and LSMS
  char physicalSystemCommandLine[256]      {};

  // Monte Carlo move set for QE systems
  int  QEMCMoveSet           {-1};

  // These are physical system specific:
  unsigned int spinModelLatticeSize  {1};                    // setting default to 1
  unsigned int spinModelDimension    {1};
  unsigned int spinConfigInitMethod  {0};
  int          numAtoms              {-1};                   // code should exit if not specified in input file (Jun 24, 17)

};

extern SimulationInfo simInfo;

#endif
