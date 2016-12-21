#ifndef GLOBALS_HPP
#define GLOBALS_HPP


struct SimulationInfo {

  int  restartFlag     {0};
  int  algorithm       {2};
  int  system          {2};
  int  rngSeed         {-1};
  char inputFile[100]  {};

  // these are system specific:
  int  size            {3};

};


#endif
