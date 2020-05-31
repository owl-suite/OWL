#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "CrystalStructure3D.hpp"
#include "Utilities/RandomNumberGenerator.hpp"


CrystalStructure3D::CrystalStructure3D(const char* inputFile, const char* spinConfigFile, int initial) : lattice(inputFile)
{

  printf("Simulation for customized 3D crystal structure: %dx%dx%d unit cells \n", lattice.unitCellDimensions[0], lattice.unitCellDimensions[1], lattice.unitCellDimensions[2]);

  assert (lattice.totalNumberOfAtoms > 0);
  spin.resize(lattice.totalNumberOfAtoms);

  if (spinConfigFile != NULL)
    readSpinConfigFile(spinConfigFile);
  else
    initializeSpinConfiguration(initial);

  // Initialize nearest neighbor lists for each atom
  neighborList.resize(lattice.totalNumberOfAtoms);
  constructPrimaryNeighborList();
  mapPrimaryToAllNeighborLists();

  //for (unsigned int i=0; i<lattice.totalNumberOfAtoms; i++)
  //  neighborList[i] = constructNeighborListFromNeighboringUnitCells(i);



  // ------OK up to here

  initializeObservables(4);      // observables[0] : total energy
                                 // observables[1] : magnetization in x-direction
                                 // observables[2] : magnetization in y-direction
                                 // observables[3] : magnetization in z-direction
  firstTimeGetMeasures = true;
  getObservables();

}


// OK
CrystalStructure3D::~CrystalStructure3D()
{

  //delete spin;
  deleteObservables();

  printf("CrystalStructure3D finished\n");

}


//void CrystalStructure3D::readCommandLineOptions()
//{ };


// OK
void CrystalStructure3D::writeConfiguration(int format, const char* filename)
{

  FILE* f;
  if (filename != NULL) f = fopen(filename, "w");
  else f = stdout;

  switch (format) {

  default : {

    fprintf(f, "\n");
    fprintf(f, "Customized 3D crystal structure: %dx%dx%d unit cells \n", lattice.unitCellDimensions[0], lattice.unitCellDimensions[1], lattice.unitCellDimensions[2]);
    fprintf(f, "Total number of atoms: %u \n", lattice.totalNumberOfAtoms);
    fprintf(f, "Measures:");
    for (unsigned int i = 0; i < numObservables; i++)
      fprintf(f, " %10.5f", observables[i]);
    fprintf(f, "\n");

    for (unsigned int i = 0; i < lattice.totalNumberOfAtoms; i++)
          fprintf(f, "%8.5f %8.5f %8.5f\n", spin[i].x, spin[i].y, spin[i].z);

  }

  }

  if (filename != NULL) fclose(f);

}


// To implement
void CrystalStructure3D::getObservablesFromScratch() 
{

/*
  //printf("!!! CALLING getObservablesFromScratch !!! \n");

  int i, j, k;
  int xLeft, yBelow, zBackward;

  // Uncomment this when observables[] are used
  //resetObservables();

  ObservableType tempE = 0.0;
  ObservableType tempMx = 0.0;
  ObservableType tempMy = 0.0;
  ObservableType tempMz = 0.0;

  for (i = 0; i < Size; i++) {
    if (i != 0) xLeft = i - 1; else xLeft = Size - 1;
    for (j= 0; j < Size; j++) {
      if (j != 0) yBelow = j - 1; else yBelow = Size - 1;
      for (k = 0; k < Size; k++) {
        if (k != 0) zBackward = k - 1; else zBackward = Size - 1;
        //observables[0] += spin[x][y].x * (spin[xLeft][y].x + spin[x][yBelow].x) + 
        //               spin[x][y].y * (spin[xLeft][y].y + spin[x][yBelow].y) +
        //               spin[x][y].z * (spin[xLeft][y].z + spin[x][yBelow].z);
        //observables[1] += spin[x][y].x;
        //observables[2] += spin[x][y].y;
        //observables[3] += spin[x][y].z;
        tempE  += spin[i][j][k].x * (spin[xLeft][j][k].x + spin[i][yBelow][k].x + spin[i][j][zBackward].x) + 
                  spin[i][j][k].y * (spin[xLeft][j][k].y + spin[i][yBelow][k].y + spin[i][j][zBackward].y) +
                  spin[i][j][k].z * (spin[xLeft][j][k].z + spin[i][yBelow][k].z + spin[i][j][zBackward].z);
        tempMx += spin[i][j][k].x;
        tempMy += spin[i][j][k].y;
        tempMz += spin[i][j][k].z;
      }
    }
  }
  //observables[0] = -observables[0];   // ferromagnetic (FO) coupling
  tempE = -tempE;

  if ((std::abs(tempE) - std::abs(observables[0])) > 10e-8) printf("Problem! tempE - observables[0] = %15.10f\n", tempE-observables[0]);
  if ((std::abs(tempMx) - std::abs(observables[1])) > 10e-8) printf("Problem! tempMx - observables[1] = %15.10f\n", tempMx-observables[1]);
  if ((std::abs(tempMy) - std::abs(observables[2])) > 10e-8) printf("Problem! tempMy - observables[2] = %15.10f\n", tempMy-observables[2]);
  if ((std::abs(tempMz) - std::abs(observables[3])) > 10e-8) printf("Problem! tempMz - observables[3] = %15.10f\n", tempMz-observables[3]);

*/
}


// To implement
void CrystalStructure3D::getObservables() 
{


}

// OK
void CrystalStructure3D::doMCMove()
{

  // Need this here since resetObservables() is not called if firstTimeGetMeasures = false
  // for (int i = 0; i < numObservables; i++)
  //   oldObservables[i] = observables[i];

  currentPosition = getUnsignedIntRandomNumber() % lattice.totalNumberOfAtoms;
  currentSpin = spin[currentPosition];

  assignRandomSpinConfiguration(currentPosition);

}


/*
void CrystalStructure3D::undoMCMove()
{
  spin[CurX][CurY][CurZ] = currentSpin;
  restoreObservables();
}
*/

// OK
void CrystalStructure3D::acceptMCMove()
{
  // update "old" observables
  for (int i=0; i<numObservables; i++)
    oldObservables[i] = observables[i];
}

// OK
void CrystalStructure3D::rejectMCMove()
{
  spin[currentPosition] = currentSpin;
  for (int i=0; i<numObservables; i++)
    observables[i] = oldObservables[i];
}

/*
void CrystalStructure3D::buildMPIConfigurationType()
{
}
*/


// To implement
void CrystalStructure3D::readSpinConfigFile(const char* spinConfigFile)
{

}



// OK
void CrystalStructure3D::initializeSpinConfiguration(int initial)
{

  for (unsigned int atomID = 0; atomID < lattice.totalNumberOfAtoms; atomID++) {

    switch (initial) {
      case 1 : {
        spin[atomID].x = 1.0;
        spin[atomID].y = 0.0;
        spin[atomID].z = 0.0;
        break;
      }
      case 2  : {
        spin[atomID].x = 0.0;
        spin[atomID].y = 1.0;
        spin[atomID].z = 0.0;
	    break;
      }
      case 3  : {
        spin[atomID].x = 0.0;
        spin[atomID].y = 0.0;
        spin[atomID].z = 1.0;
	      break;
      }
      case 4  : {
        spin[atomID].x = 0.0;
        spin[atomID].y = 0.0;
        spin[atomID].z = (atomID % 2 == 0) ? 1.0 : -1.0;
        break;
      }
      default  : {
        assignRandomSpinConfiguration(atomID);
      }
    }

  }

}


// OK
void CrystalStructure3D::assignRandomSpinConfiguration(unsigned int currentAtom)
{

  double r1, r2, rr;

  do {
    r1 = 2.0 * getRandomNumber();
    r2 = 2.0 * getRandomNumber();                
    rr = r1 * r1 + r2 * r2;
  } while (rr > 1.0);

  spin[currentAtom].x = 2.0 * r1 * sqrt(1.0 - rr);
  spin[currentAtom].y = 2.0 * r2 * sqrt(1.0 - rr);
  spin[currentAtom].z = 1.0 - 2.0 * rr;

}


// To implement
void CrystalStructure3D::readHamiltonianTerms(const char* inputFile)
{

  
}


// To implement
ObservableType CrystalStructure3D::nearestNeighborInterations()
{

  return 0.0;

}


// Construct the neighbor list for each atom, given the list of neighboring unit cells to check
std::vector<NeighboringAtomInfo> CrystalStructure3D::constructNeighborListFromNeighboringUnitCells(unsigned int atom1_global)
{

    // A neighbor list where pairwise interactions with current atom will be checked
    std::vector<NeighboringAtomInfo> atomList;

    // Current unit cell is always (0,0,0) relatively.
    unsigned int localUnitCellIndex = lattice.getRelativeUnitCellIndex(0, 0, 0);
    unsigned int globalUnitCellIndex =  atom1_global / lattice.unitCell.number_of_atoms;
    
    unsigned int atom1 = lattice.getAtomIndex(localUnitCellIndex,  atom1_global % lattice.unitCell.number_of_atoms);

    for (unsigned int j=0; j<lattice.nearestNeighborUnitCellList[globalUnitCellIndex].size(); j++) {
      for (unsigned int k=0; k<lattice.unitCell.number_of_atoms; k++) {

        unsigned int atom2 = lattice.getAtomIndex(j, k);
        unsigned int atom2_global = lattice.getAtomIndex(lattice.nearestNeighborUnitCellList[globalUnitCellIndex][j], k);
        if (atom1_global == atom2_global) break;
        double distance = lattice.getRelativePairwiseDistance(atom1, atom2);

        if (distance <= 1.0)
          atomList.push_back({atom2_global, distance, 0.0, 0.0});
        
      }
    }

    std::sort(atomList.begin(), atomList.end(), 
              [](const auto& a, const auto& b) { return a.atomID < b.atomID; }
    );

    std::cout << "Sorted neighbor list of " << atom1_global << ":\n";
    std::cout << "Atom    distance \n";
    for (auto i : atomList)
      std::cout << i.atomID << " " << i.distance << "\n";

  return atomList;

}


// Construct the neighbor list for each atom in a unit cell
// Unit cell and atomID are relative to the reference unit cell
void CrystalStructure3D::constructPrimaryNeighborList()
{

  double distance {0.0};
  double J_ij     {0.0};          // Exchange coupling
  double D_ij     {0.0};          // Dzyaloshinskii-Moriya (DM) interaction
  double dx       {0.0};
  double dy       {0.0};
  double dz       {0.0};

  primaryNeighborList.resize(lattice.unitCell.number_of_atoms);
  unsigned int localUnitCellIndex = lattice.getRelativeUnitCellIndex(0, 0, 0);

  for (unsigned int atomID=0; atomID<lattice.unitCell.number_of_atoms; atomID++) {
    unsigned int atom1 = lattice.getAtomIndex(localUnitCellIndex,  atomID);

    for (unsigned int j=0; j<lattice.nearestNeighborUnitCellList[localUnitCellIndex].size(); j++) {
      for (unsigned int k=0; k<lattice.unitCell.number_of_atoms; k++) {
        unsigned int atom2 = lattice.getAtomIndex(j, k);
        if (atom1 == atom2) break;
        dx = lattice.relativeAtomicPositions(0, atom2) - lattice.relativeAtomicPositions(0, atom1);
        dy = lattice.relativeAtomicPositions(1, atom2) - lattice.relativeAtomicPositions(1, atom1);
        dz = lattice.relativeAtomicPositions(2, atom2) - lattice.relativeAtomicPositions(2, atom1);
        distance = lattice.getRelativePairwiseDistance(atom1, atom2);

        std::cout << "atom " << atom1 << " (" << lattice.relativeAtomicPositions(0, atom1) << ", "
                                              << lattice.relativeAtomicPositions(1, atom1) << ", "
                                              << lattice.relativeAtomicPositions(2, atom1) 
                  << ") , atom " << atom2 << " (" << lattice.relativeAtomicPositions(0, atom2) << ", "
                                                  << lattice.relativeAtomicPositions(1, atom2) << ", "
                                                  << lattice.relativeAtomicPositions(2, atom2) 
                  << ") from unit cell (" 
                  << lattice.relativeUnitCellVectors(0,j) << " " 
                  << lattice.relativeUnitCellVectors(1,j) << " " 
                  << lattice.relativeUnitCellVectors(2,j) << ") : " 
                  << dx << " , " << dy << " , " << dz << " . " << distance << "\n";

        J_ij = assignExchangeCouplings(dx, dy, dz, distance);
        D_ij = assignDMInteractions(dx, dy, dz, distance);

        if (distance <= 1.0)
          primaryNeighborList[atomID].push_back({atom2, distance, J_ij, D_ij});
      }
    }
    
    std::cout << "Primary neighbor list of " << atom1 << ":\n";
    std::cout << "Atom    distance \n";
    for (auto i : primaryNeighborList[atomID])
      std::cout << i.atomID << " " << i.distance << " " << i.J_ij << " " << i.D_ij << "\n";

  }

  std::cout << "CrystalStructure3D: Constructed primary neighbor lists for all atoms in a unit cell. \n";

}


void CrystalStructure3D::mapPrimaryToAllNeighborLists()
{

  unsigned int thisAtom, atom_tmp, relative_uc, real_uc, atomID_in_uc, atomID;

  for (unsigned int i=0; i<lattice.numberOfUnitCells; i++) {
    for (unsigned int j=0; j<lattice.unitCell.number_of_atoms; j++) {

      thisAtom = lattice.getAtomIndex(i,j);

      for (auto k : primaryNeighborList[j]) {
        atom_tmp = k.atomID;

        relative_uc = atom_tmp / lattice.unitCell.number_of_atoms;       // which relative unit cell the neighoring atom in?
        atomID_in_uc = atom_tmp % lattice.unitCell.number_of_atoms;      // which atom is the neighoring atom in a unit cell?
        real_uc = lattice.nearestNeighborUnitCellList[i][relative_uc];
        atomID = lattice.getAtomIndex(real_uc, atomID_in_uc);

        neighborList[thisAtom].push_back({atomID, k.distance, k.J_ij, k.D_ij});
        //std::cout << atomID << " " << k.distance << "\n";
      }

      std::sort(neighborList[thisAtom].begin(), neighborList[thisAtom].end(), 
                [](const auto& a, const auto& b) { return a.atomID < b.atomID; }
      );

      std::cout << "Mapped Neighbor list of " << thisAtom << ":\n";
      std::cout << "Atom     distance      J_ij       D_ij \n";
      for (auto m : neighborList[thisAtom] )
        std::cout << m.atomID << " " << m.distance << " " << m.J_ij << " " << m.D_ij << "\n";

    }
  }

  std::cout << "CrystalStructure3D: Mapped primary neighbor lists to all atoms. \n";

}

// Note: Coupling measured in meV
double CrystalStructure3D::assignExchangeCouplings(double dx, double dy, double dz, double dr)
{

  const double dr_ref1   {0.613054};     // should be defined in input file if possible
  const double dr_ref2   {0.913759};
  const double dr_ref3   {0.957453};
  const double dr_ref4   {1.0};
  const double ref1      {0.22956};
  const double ref2      {0.27044};
  const double ref3      {0.5};
  const double threshold {0.0001};
  double coupling        {0.0};

  auto sameMagitude = [=](double a, double b) -> bool {return fabs(fabs(a) - fabs(b)) < threshold; }; 

  if (sameMagitude(dr, dr_ref1)) {         // nearest-neighbor coupling
    if (sameMagitude(dx, ref1))             coupling = 5.746999335;
    else if (sameMagitude(dx, ref2))        coupling = 5.746992166;
    else if (sameMagitude(dx, ref3))        coupling = 5.750420632;
  }  
  else if (sameMagitude(dr, dr_ref2)) {    // next nearest-neighbor coupling
    if (sameMagitude(dx, ref1))             coupling = -1.288546897;
    else if (sameMagitude(dx, ref3))        coupling = -1.286718254;
    else if (sameMagitude(dx, ref1+ref3))   coupling = -1.285652206;
  }  
  else if (sameMagitude(dr, dr_ref3)) {    // next next nearest-neighbor coupling
    if (sameMagitude(dx, ref2))             coupling = -1.019740021;
    else if (sameMagitude(dx, ref3))        coupling = -1.021497653;
    else if (sameMagitude(dx, ref2+ref3))   coupling = -1.019455211;
  }  
  else if (sameMagitude(dr, dr_ref4)) {    // 4th nearest-neighbor coupling
    if (sameMagitude(dx, dr_ref4))          coupling = 0.231944895;
    else if (sameMagitude(dy, dr_ref4))     coupling = 0.233575161;
    else if (sameMagitude(dz, dr_ref4))     coupling = 0.234596638;
  }

  return coupling;

}


double CrystalStructure3D::assignDMInteractions(double dx, double dy, double dz, double dr)
{
  
  const double dr_ref1   {0.613054};     // should be defined in input file if possible
  const double dr_ref2   {0.913759};
  const double dr_ref3   {0.957453};
  const double dr_ref4   {1.0};
  const double ref1      {0.22956};
  const double ref2      {0.27044};
  const double ref3      {0.5};
  const double threshold {0.001};
  double coupling        {0.0};
  
  auto sameMagitude = [=](double a, double b) -> bool { return fabs(fabs(a) - fabs(b)) < threshold; };
  auto sameSign     = [=](double a, double b) -> bool { return a * b > 0.0; };

  if (sameMagitude(dr, dr_ref1)) {         // nearest-neighbor coupling
    if (sameMagitude(dx, ref1))             coupling = 0.17494;
    else if (sameMagitude(dx, ref3))        coupling = -0.06926;
    else if (sameMagitude(dx, ref2)) {
      if (sameSign(dx, dz))                 coupling = -0.09617;
      else                                  coupling = 0.09617;
    } 
  }
  else if (sameMagitude(dr, dr_ref2)) {    // next nearest-neighbor coupling
    if (sameMagitude(dx, ref1))             coupling = 0.05245;
    else if (sameMagitude(dx, ref3))        coupling = -0.04845;
    else if (sameMagitude(dx, ref1+ref3)) {
      if (sameSign(dx, dz))                 coupling = -0.06946;
      else                                  coupling = 0.06946;
    }
  }
  else if (sameMagitude(dr, dr_ref3)) {    // next next nearest-neighbor coupling
    if (sameMagitude(dx, ref3))             coupling = -0.02081;
    else if (sameMagitude(dx, ref2+ref3))   coupling = 0.02138;
    else if (sameMagitude(dx, ref2)) {
      if (sameSign(dx, dz))                 coupling = 0.09933;
      else                                  coupling = -0.09933;
    }
  }
  else if (sameMagitude(dr, dr_ref4)) {    // 4th nearest-neighbor coupling
    if (sameMagitude(dx, dr_ref4))          coupling = 1.0;
    else if (sameMagitude(dy, dr_ref4))     coupling = 1.0;
    else if (sameMagitude(dz, dr_ref4))     coupling = 1.0;
  }

  return coupling;
  
}