#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Matrix.hpp"


// Moved to RandomNumberGenerator.hpp
//void initializeRandomNumberGenerator(int seed = -1);
//inline double getRandomNumber()
//{
//  return (double(rand())/double(RAND_MAX) - 0.5);
//}



// "Private"
void writeAtomicPositions(Matrix<double> atom_positions);

void writeLatticeVectors(Matrix<double> cell_vectors);

// "Public"
void proposeMCmoves(Matrix<double> &atom_positions, Matrix<double> &cell_vectors);


