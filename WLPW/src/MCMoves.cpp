#include <cmath>
#include <limits>
#include "MCMoves.hpp"


void initializeRandomNumberGenerator(int seed)
{
  /* initialize random seed */
  if (seed == -1) {
    srand(time(NULL));
    std::cout << "No random number seed supplied. Take current time as a seed." << std::endl;
  }
  else {
    std::cout << "Random number seed supplied: " << seed << std::endl;
    srand(seed);
  }
  
  //std::cout << "RAND_MAX = " << RAND_MAX << std::endl;
}


void writeAtomicPositions(Matrix<double> atom_positions)
{

  printf("YingWai's debug: atom_position inside displaceAnAtom\n");
  for(int i=0; i<atom_positions.n_col(); i++) {
    printf("atom %d : ", i);
    for(int j=0; j<atom_positions.n_row(); j++)
      printf(" %14.9f ", atom_positions(j,i));
    printf("\n");
  }

}


void writeLatticeVectors(Matrix<double> cell_vectors)
{

  printf("YingWai's debug: lattice_vector inside stretchCrystalCell\n");
  for(int i=0; i<cell_vectors.n_col(); i++) {
    printf("vector %d : ", i);
    for(int j=0; j<cell_vectors.n_row(); j++)
      printf(" %14.9f ", cell_vectors(j,i));
    printf("\n");
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
    dr = 2.0 * getRandomNumber() * dr_max;
    // new atomic positions
    atom_positions(i,ranAtom) += dr;
  }

  // YingWai's debugging check
  //writeAtomicPositions(atom_positions);

}



void stretchCrystalCell(Matrix<double> &cell_vectors)
{

  // YingWai's note (data layout of cell_vectors[i][j]):
  //    j -->
  //  i    a_x b_x c_x
  //  |    a_y b_y c_y
  //  v    a_z b_z c_z

  // Define maximum stretching magnitude in Angstrom  (should be moved to .h later)
  double dL_max = 0.1;

  // Choose a vector randomly to stretch
  int ranLatticeVector = rand() % cell_vectors.n_row() ;

  // Lengths of lattice vector
  double L = 0.0;
  for (int j=0; j<cell_vectors.n_col(); j++)
    L += cell_vectors(j,ranLatticeVector) * cell_vectors(j,ranLatticeVector);
  L = sqrt(L);

  // Angles of the lattice vectors with the Cartesian coordinate system (x,y,z)
  // x = L * sin(theta) * cos(phi)
  // y = L * sin(theta) * sin(phi)
  // z = L * cos(theta)
  double cos_theta = cell_vectors(2,ranLatticeVector) / L;
  double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
  double cos_phi, sin_phi;

  if (sin_theta < std::numeric_limits<double>::epsilon()) {   // the vector is in z-direction
    cos_phi = 0.0;
    sin_phi = 0.0;
  }
  else {
    cos_phi  = cell_vectors(0,ranLatticeVector) / (L * sin_theta);
    sin_phi  = cell_vectors(1,ranLatticeVector) / (L * sin_theta);
  }

  // Randomly choose a stretching magnitude
  double dL = 2.0 * getRandomNumber() * dL_max;
  L += dL;

  // Recalculate cell vector's components with new length
  cell_vectors(0,ranLatticeVector) = L * sin_theta * cos_phi;
  cell_vectors(1,ranLatticeVector) = L * sin_theta * sin_phi;
  cell_vectors(2,ranLatticeVector) = L * cos_theta;

  //YingWai's debugging check
  //writeLatticeVectors(cell_vectors);

}


void proposeMCmoves(Matrix<double> &atom_positions, Matrix<double> &cell_vectors)
{

  enum MCMoves {displace_atom, stretch_crystal_cell};
  //std::cout << "enum DisplaceAtom = " << displace_atom << std::endl;
  //std::cout << "enum ElongateCrystalCell = " << stretch_crystal_cell << std::endl;

  std::cout << "size = " << atom_positions.size()  << std::endl;
  std::cout << "nRow = " << atom_positions.n_row() << std::endl;
  std::cout << "nCol = " << atom_positions.n_col() << std::endl;
  std::cout << "lDim = " << atom_positions.l_dim() << std::endl;
  std::cout << " i    atom_positions " << std::endl;

  //Choose a MC move to perform
  int r = rand() % 2;
  switch(r)
  {
    case 0 :
      displaceAnAtom(atom_positions);
      std::cout << "MCMove: displace an atom\n";
      break;
    case 1 :
      stretchCrystalCell(cell_vectors);
      std::cout << "MCMove: stretch crystal cell\n";
      break;
    default:
      std::cout << "MCMove: invalid choice. r = %d " << r << std::endl;
      break;
  }
  
}
