#include "MCAlgorithms.hpp"

void YingWaisCheck()
{

  initializeRandomNumberGenerator();
  Matrix<double> pos_array_trial;   // Set up the position array (in angstrom)
  pos_array_trial.resize(5,10);
  Matrix<double> cell_array_trial;  // Set up the cell vector array (in angstrom)
  cell_array_trial.resize(3,3);
  for(size_t i=0; i<50; i++)
    pos_array_trial[i] = i;
  for(size_t i=0; i<9; i++)
    cell_array_trial[i] = i;
  proposeMCmoves(pos_array_trial, cell_array_trial);
  
  //for(size_t i=0; i<50; i++)
  //  cout << "main: " << i << "  " << pos_array_trial[i] << endl;

}
