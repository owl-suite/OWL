#include <iostream>
#include <cstdlib>
#include <ctime>

#include "Matrix.hpp"
using namespace std;


void initializeRandomNumberGenerator(int seed = -1);

inline double getRandomNumber()
{
  return (rand()/RAND_MAX - 0.5);
}

void proposeMCmoves(Matrix<double> &atom_positions, Matrix<double> &cell_vectors);
