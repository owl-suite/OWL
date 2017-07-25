#include <cstdio>
#include <cmath>
#include "MCAlgorithms.hpp"


MonteCarloAlgorithm::MonteCarloAlgorithm()
{
  printf("MonteCarloAlgorithm constructor:\n");
  //totalMCsteps  = 0;
  //acceptedMoves = 0;
  //rejectedMoves = 0;
  printf("--totalMCsteps:  %lu \n", totalMCsteps);
  printf("--acceptedMoves: %lu \n", acceptedMoves);
  printf("--rejectedMoves: %lu \n", rejectedMoves);
}


/*
MonteCarloAlgorithm::~MonteCarloAlgorithm()
{

}
*/

