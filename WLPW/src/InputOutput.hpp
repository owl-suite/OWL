# include <iostream>
# include <fstream>
# include <cstdio>
# include "Matrix.hpp"

void writeEnergyFile(char fileName[], double energy);

void writeQErestartFile(char fileName[], Matrix<double> atom_postiion, Matrix<double> cell_vectors);
