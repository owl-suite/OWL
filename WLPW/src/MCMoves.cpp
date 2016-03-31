#include <cmath>
#include "MCMoves.hpp"


void initializeRandomNumberGenerator(int seed)
{
  /* initialize random seed */
  if (seed == -1) {
    srand(time(NULL));
    cout << "No random number seed supplied. Take current time as a seed." << endl;
  }
  else {
    cout << "Random number seed supplied: " << seed << endl;
    srand(seed);
  }
}


// Takes in old atomic positions and modify
void displaceAnAtom(Matrix<double> &atom_positions)
{
  // YingWai's note (data layout of Matrix class):
  // (i) mat(i,j)
  //     resize(nRow,nCol); nRow = dim. of coordinates, nCol = # of atoms
  //       j ->
  //     i    x0  x1 ... x[j-1]
  //     |    y0  y1     y[j-1]
  //     v    z0  z1 ... z[j-1]
  //
  // (ii) mat[k]
  //      index is calculated by k = j*lDim + i
  //       0  1  2  3  4  5 ... 
  //      x0 y0 z0 x1 y1 z1 ... 


  // Define maximum displacement magnitude in Angstrom  (should be moved to .h later)
  double dr_max = 0.1;

  // Displacement
  double dr = 0.0;

  // Choose an atom randomly
  int ranAtom = rand() % atom_positions.n_col() ;

  // For each of the x-,y-,z-direction,
  for (int i=0; i<atom_positions.n_row(); i++)
  {
    // randomly choose a displacement magnitude
    dr = 2.0 * (rand()/RAND_MAX - 0.5) * dr_max;
    // new atomic positions
    atom_positions(i,ranAtom) += dr;
  }

}



void stretchCrystalCell(Matrix<double> &cell_vectors)
{

  // YingWai's note (data layout of cell_vectors[i][j]):
  //    j -->
  //  i    a_x a_y a_z
  //  |    b_x b_y b_z
  //  v    c_x c_y c_z

  // Define maximum stretching magnitude in Angstrom  (should be moved to .h later)
  double dL_max = 0.1;

  // Choose a vector randomly to stretch
  int ranLatticeVector = rand() % cell_vectors.n_row() ;

  // Lengths of lattice vector
  double L = 0.0;
  for (int j=0; j<cell_vectors.n_col(); j++)
    L += cell_vectors(ranLatticeVector,j) * cell_vectors(ranLatticeVector,j);
  L = sqrt(L);

  // Angles of the lattice vectors with the Cartesian coordinate system (x,y,z)
  // x = L * sin(theta) * cos(phi)
  // y = L * sin(theta) * sin(phi)
  // z = L * cos(theta)
  double cos_theta = cell_vectors(ranLatticeVector,2) / L;
  double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
  double cos_phi   = cell_vectors(ranLatticeVector,0) / (L * sin_theta);
  double sin_phi   = cell_vectors(ranLatticeVector,1) / (L * sin_theta);

  // Randomly choose a stretching magnitude
  double dL = 2.0 * (rand()/RAND_MAX - 0.5) * dL_max;
  L += dL;

  // Recalculate cell vector's components with new length
  cell_vectors(ranLatticeVector,0) = L * sin_theta * cos_phi;
  cell_vectors(ranLatticeVector,1) = L * sin_theta * sin_phi;
  cell_vectors(ranLatticeVector,2) = L * cos_theta;

}



void proposeMCmoves(Matrix<double> &atom_positions, Matrix<double> &cell_vectors)
{

  enum MCMoves {displace_atom, stretch_crystal_cell};
  cout << "enum DisplaceAtom = " << displace_atom << endl;
  cout << "enum ElongateCrystalCell = " << stretch_crystal_cell << endl;

  cout << "size = " << atom_positions.size()  << endl;
  cout << "nRow = " << atom_positions.n_row() << endl;
  cout << "nCol = " << atom_positions.n_col() << endl;
  cout << "lDim = " << atom_positions.l_dim() << endl;
  cout << " i    atom_positions " << endl;

/*
  for (size_t i=0; i<atom_positions.size(); i++) {
    cout << i << "  " << atom_positions[i] << endl;
    atom_positions[i] += 1.1;
  }
*/

  //Choose a MC move to perform
  int r = rand() % 2;
  switch(r)
  {
    case 0 :
      displaceAnAtom(atom_positions);
      cout << "MCMove: displace an atom\n";
      break;
    case 1 :
      stretchCrystalCell(cell_vectors);
      cout << "MCMove: stretch crystal cell\n";
      break;
    default:
      cout << "MCMove: invalid choice. r = %d " << r << endl;
      break;
  }
  
}
