# include <iostream>
# include <fstream>
# include <cstdio>
# include "Matrix.hpp"

void writeSystemFile(char fileName[], double energy, Matrix<double> atom_positions, Matrix<double> cell_vectors);

void writeQErestartFile(char fileName[], Matrix<double> atom_positions, Matrix<double> cell_vectors);
