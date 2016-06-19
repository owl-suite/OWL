#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Matrix.hpp"


// "Private"
void initializeRandomNumberGenerator(int seed = -1);

void writeAtomicPositions(Matrix<double> atom_positions);

void writeLatticeVectors(Matrix<double> cell_vectors);

inline double getRandomNumber()
{
  return (double(rand())/double(RAND_MAX) - 0.5);
}

// "Public"
void proposeMCmoves(Matrix<double> &atom_positions, Matrix<double> &cell_vectors);


