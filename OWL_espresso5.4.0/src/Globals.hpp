//  This head file stores all the global variables / data structures
//  for sharing among different part of the code

#ifndef GLOBALS_HPP
#define GLOBALS_HPP


struct SimulationInfo
{

  // Basic simulation info
  int  restartFlag           {0};
  int  algorithm             {2};
  int  system                {2};
  int  rngSeed               {-1};
  char MCInputFile[255]      {};

  // MPI info
  int  numWalkers            {1};
  int  numMPIranksPerWalker  {1};

  // Command line options for QE and LSMS
  char physicalSystemCommandLine[256]      {};

  // these are system specific:
  int  spinModelLatticeSize  {3};   // considering setting default to -1...
  int  numAtoms              {3};   // it should exit if not specified in input file (Jun 24, 17)

};


#endif
