#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <filesystem>
#include <fstream>
#include "CrystalStructure3D.hpp"
#include "Utilities/CheckFile.hpp"
#include "Utilities/RandomNumberGenerator.hpp"


CrystalStructure3D::CrystalStructure3D(const char* inputFile, int initial) : lattice(inputFile)
{

  printf("\n");
  printf("Simulation for customized 3D crystal structure: %dx%dx%d unit cells \n", 
         lattice.unitCellDimensions[0], lattice.unitCellDimensions[1], lattice.unitCellDimensions[2]);

  assert (lattice.totalNumberOfAtoms > 0);
  setSystemSize(lattice.totalNumberOfAtoms);
  spin.resize(systemSize);
  
  // TODO: This should be incorporated into the constructor of the Hamiltonian class later when it is implemented,
  //       together with the reading of Hamiltonian terms. (July 7, 20)
  if (std::filesystem::exists(inputFile))
    readInteractionCutoffDistance(inputFile);
  //else  (TODO)
  //  printf("   Input file '%s' not found. Interaction cutoff distance will be set to nearest-neighbor only.\n", inputFile);
  
  // Initialize nearest neighbor lists for each atom
  neighborList.resize(systemSize);
  constructPrimaryNeighborList();
  mapPrimaryToAllNeighborLists();

  //for (unsigned int i=0; i<systemSize; i++)
  //  neighborList[i] = constructNeighborListFromNeighboringUnitCells(i);
  
  initializeObservables(6); 
  observableName.push_back("Total energy, E");                            // observables[0] : total energy
  observableName.push_back("Magnetization in x-direction, M_x");          // observables[1] : magnetization in x-direction
  observableName.push_back("Magnetization in y-direction, M_y");          // observables[2] : magnetization in y-direction
  observableName.push_back("Magnetization in z-direction, M_z");          // observables[3] : magnetization in z-direction
  observableName.push_back("Total magnetization, M");                     // observables[4] : total magnetization
  observableName.push_back("4th order magnetization, M^4");               // observables[5] : total magnetization to the order 4

  // Initialize configuration from file if applicable
  if (std::filesystem::exists("config_initial.dat"))
    readSpinConfigFile("config_initial.dat");
  else if (simInfo.restartFlag && std::filesystem::exists("configurations/config_checkpoint.dat"))
    readSpinConfigFile("configurations/config_checkpoint.dat");
  else
    initializeSpinConfiguration(initial);

  firstTimeGetMeasures = true;
  getObservablesFromScratch();

}


// OK
CrystalStructure3D::~CrystalStructure3D()
{

  //delete spin;
  deleteObservables();

  printf("Exiting CrystalStructure3D class...\n");

}


//void CrystalStructure3D::readCommandLineOptions()
//{ };



void CrystalStructure3D::writeConfiguration(int format, const char* filename)
{

  FILE* configFile;
  if (filename != NULL) configFile = fopen(filename, "w");
  else configFile = stdout;

  switch (format) {

  default : {

    fprintf(configFile, "# Customized 3D crystal structure: %dx%dx%d unit cells \n\n", lattice.unitCellDimensions[0], lattice.unitCellDimensions[1], lattice.unitCellDimensions[2]);
    fprintf(configFile, "TotalNumberOfAtoms %u \n", systemSize);
    fprintf(configFile, "Observables ");
    for (unsigned int i = 0; i < numObservables; i++)
      fprintf(configFile, " %15.8f", observables[i]);
    fprintf(configFile, "\n\n");

    fprintf(configFile, "SpinConfiguration\n");
    for (unsigned int i = 0; i < systemSize; i++)
      fprintf(configFile, "%8.5f %8.5f %8.5f\n", spin[i].x, spin[i].y, spin[i].z);

  }

  }

  if (filename != NULL) fclose(configFile);

}


void CrystalStructure3D::getObservablesFromScratch() 
{

  observables[0] = getExchangeInterations();
  //observables[0] = getExchangeInterations() + getDzyaloshinskiiMoriyaInterations();
  std::tie(observables[1], observables[2], observables[3], observables[4]) = getMagnetization();
  observables[5] = pow(observables[4], 4.0);

  firstTimeGetMeasures = false;

}


void CrystalStructure3D::getObservables() 
{

  observables[0] += getDifferenceInExchangeInterations();
  //observables[0] += getDifferenceInExchangeInterations() + getDifferenceInDzyaloshinskiiMoriyaInterations();
  observables[1] += spin[currentPosition].x - oldSpin.x;
  observables[2] += spin[currentPosition].y - oldSpin.y;
  observables[3] += spin[currentPosition].z - oldSpin.z;
  ObservableType temp = observables[1] * observables[1] + observables[2] * observables[2] + observables[3] * observables[3];
  observables[4] = sqrt(temp);
  observables[5] = temp * temp;

}



void CrystalStructure3D::doMCMove()
{

  // Need this here since resetObservables() is not called if firstTimeGetMeasures = false
  // for (unsigned int i = 0; i < numObservables; i++)
  //   oldObservables[i] = observables[i];

  currentPosition = getUnsignedIntRandomNumber() % systemSize;
  oldSpin = spin[currentPosition];

  assignRandomSpinDirection(currentPosition);

}



/*
void CrystalStructure3D::undoMCMove()
{
  spin[CurX][CurY][CurZ] = oldSpin;
  restoreObservables();
}
*/



void CrystalStructure3D::acceptMCMove()
{
  // update "old" observables
  for (unsigned int i = 0; i < numObservables; i++)
    oldObservables[i] = observables[i];
}



void CrystalStructure3D::rejectMCMove()
{
  spin[currentPosition] = oldSpin;
  for (unsigned int i = 0; i < numObservables; i++)
    observables[i] = oldObservables[i];
}

/*
void CrystalStructure3D::buildMPIConfigurationType()
{
}
*/


void CrystalStructure3D::readSpinConfigFile(const std::filesystem::path& spinConfigFile)
{
  std::cout << "\n   CrystalStructure3D class reading configuration file: " << spinConfigFile << "\n";

  std::ifstream inputFile(spinConfigFile);
  std::string line, key;
  unsigned int numberOfAtoms {0};

  if (inputFile.is_open()) {

    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
          std::istringstream lineStream(line);
          lineStream >> key;

          if (key.compare(0, 1, "#") != 0) {

            if (key == "TotalNumberOfAtoms") {
              lineStream >> numberOfAtoms;
              //std::cout << "   CrystalStructure3D: numberOfAtoms = " << numberOfAtoms << "\n";
              continue;
            }
            else if (key == "Observables") {
              unsigned int counter = 0;
              while (lineStream && counter < numObservables) {
                lineStream >> observables[counter];
                //std::cout << "   CrystalStructure3D: observables[" << counter << "] = " << observables[counter] << "\n";
                counter++;
              }
              continue;
            }
            else if (key == "SpinConfiguration") {
              //std::cout << "   CrystalStructure3D:  Spin Configuration read: \n";
              for (unsigned int atomID=0; atomID<numberOfAtoms; atomID++) {
                lineStream.clear();
                std::getline(inputFile, line);               
                if (!line.empty()) lineStream.str(line);
                lineStream >> spin[atomID].x >> spin[atomID].y >> spin[atomID].z;
                //printf("      %8.5f %8.5f %8.5f\n", spin[atomID].x, spin[atomID].y, spin[atomID].z);
              }
              continue;
            }

          }

      }
    }

    inputFile.close();

  }

  // Sanity checks:
  assert(numberOfAtoms == systemSize);

}



void CrystalStructure3D::readInteractionCutoffDistance(const char* mainInputFile)
{
  std::cout << "   CrystalStructure3D class reading input file: " << mainInputFile << "\n";

  std::ifstream inputFile(mainInputFile);
  std::string line, key;

  if (inputFile.is_open()) {

    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
          std::istringstream lineStream(line);
          lineStream >> key;

          if (key.compare(0, 1, "#") != 0) {

            if (key == "InteractionCutoffDistance") {
              lineStream >> interactionCutoffDistance;
              std::cout << "   CrystalStructure3D: interaction cutoff distance = " << interactionCutoffDistance << "\n";
              continue;
            }

          }

      }
    }

    inputFile.close();

  }

}



void CrystalStructure3D::initializeSpinConfiguration(int initial)
{

  for (unsigned int atomID = 0; atomID < systemSize; atomID++) {

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
        assignRandomSpinDirection(atomID);
      }
    }

  }

}


void CrystalStructure3D::assignRandomSpinDirection(unsigned int currentAtom)
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


// TODO: To implement
/*
void CrystalStructure3D::readHamiltonianTerms(const char* inputFile)
{

  
}
*/


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

        if (distance <= interactionCutoffDistance)
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

  std::cout << "\n   Construct primary neighbor list:\n";

  primaryNeighborList.resize(lattice.unitCell.number_of_atoms);
  unsigned int localUnitCellIndex = lattice.getRelativeUnitCellIndex(0, 0, 0);

  for (unsigned int atomID=0; atomID<lattice.unitCell.number_of_atoms; atomID++) {
    unsigned int atom1 = lattice.getAtomIndex(localUnitCellIndex,  atomID);

    for (unsigned int j=0; j<lattice.nearestNeighborUnitCellList[localUnitCellIndex].size(); j++) {
      for (unsigned int k=0; k<lattice.unitCell.number_of_atoms; k++) {
        unsigned int atom2 = lattice.getAtomIndex(j, k);
        //if (atom1 == atom2) break;                // avoids double counting within the same unit cell
        if (atom1 == atom2) continue;               // avoids putting the reference atom itself into the neighbor list
        dx = lattice.relativeAtomicPositions(0, atom2) - lattice.relativeAtomicPositions(0, atom1);
        dy = lattice.relativeAtomicPositions(1, atom2) - lattice.relativeAtomicPositions(1, atom1);
        dz = lattice.relativeAtomicPositions(2, atom2) - lattice.relativeAtomicPositions(2, atom1);
        distance = lattice.getRelativePairwiseDistance(atom1, atom2);

        //std::cout << "atom " << atom1 << " (" << lattice.relativeAtomicPositions(0, atom1) << ", "
        //                                      << lattice.relativeAtomicPositions(1, atom1) << ", "
        //                                      << lattice.relativeAtomicPositions(2, atom1) 
        //          << ") , atom " << atom2 << " (" << lattice.relativeAtomicPositions(0, atom2) << ", "
        //                                          << lattice.relativeAtomicPositions(1, atom2) << ", "
        //                                          << lattice.relativeAtomicPositions(2, atom2) 
        //          << ") from unit cell (" 
        //          << lattice.relativeUnitCellVectors(0,j) << " " 
        //          << lattice.relativeUnitCellVectors(1,j) << " " 
        //          << lattice.relativeUnitCellVectors(2,j) << ") : " 
        //          << dx << " , " << dy << " , " << dz << " . " << distance << "\n";

        // Store the distance if it is not yet in neighborDistances
        if (!isFoundInVector(distance, neighborDistances))
          neighborDistances.push_back(distance);

        // Add the atom to neighbor list if within cutoff
        if (distance <= interactionCutoffDistance) {
          J_ij = assignExchangeCouplings(dx, dy, dz, distance);
          D_ij = assignDzyaloshinskiiMoriyaInteractions(dz, distance);
          //J_ij = assignExchangeCouplings_testing(dx, dy, dz, distance);
          primaryNeighborList[atomID].push_back({atom2, distance, J_ij, D_ij});
        }

      }
    }

    // Print the primary neighbor list for the current atom
    std::cout << "\n     Primary neighbor list of " << atom1 << ":\n";
    std::cout << "          Atom      Distance        J_ij          D_ij\n";
    for (auto i : primaryNeighborList[atomID])
      printf("     %8d   %12.6f   %12.8f   %12.8f\n", i.atomID, i.distance, i.J_ij, i.D_ij);

  }

  // Sort the neighbor distance list and print out
  std::sort(neighborDistances.begin(), neighborDistances.end(), 
            [](const auto& a, const auto& b) { return a < b; }
  );
  std::cout << "\n     Neighboring distances: \n";
  for (unsigned int i=0; i<neighborDistances.size(); i++)
    printf("     %2dth neighbor : %12.6f \n", i, neighborDistances[i]);

  std::cout << "\n   Constructed primary neighbor lists for all atoms in a unit cell. \n";

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

      //std::cout << "Mapped Neighbor list of " << thisAtom << ":\n";
      //std::cout << "Atom     distance      J_ij       D_ij \n";
      //for (auto m : neighborList[thisAtom] )
      //  std::cout << m.atomID << " " << m.distance << " " << m.J_ij << " " << m.D_ij << "\n \n";

    }
  }

  std::cout << "   Mapped primary neighbor lists to all atoms in system. \n";

}

// Note: Coupling measured in meV
// Ad hoc for our system for now.
// TODO: All the reference dr's and coupling strengths should be read from input file (using readHamiltonianTerms()).
double CrystalStructure3D::assignExchangeCouplings(double dx, double dy, double dz, double dr)
{

  const double dr_ref1   {0.613054};
  const double dr_ref2   {0.913759};
  const double dr_ref3   {0.957453};
  const double dr_ref4   {1.0};
  const double dr_ref5   {1.17296};
  const double dr_ref6   {1.35461};
  const double dr_ref7   {1.38445};

  const double ref1      {0.22956};
  const double ref2      {0.27044};
  const double ref3      {0.5};
  double coupling        {0.0};

  //auto sameMagnitude = [=](double a, double b) -> bool {return fabs(fabs(a) - fabs(b)) < threshold; }; 

  if (sameMagnitude(dr, dr_ref1)) {         // nearest-neighbor coupling
    if (sameMagnitude(dx, ref1))             coupling = 5.737246246;
    else if (sameMagnitude(dx, ref2))        coupling = 5.737112873;
    else if (sameMagnitude(dx, ref3))        coupling = 5.740528109;
  }  
  else if (sameMagnitude(dr, dr_ref2)) {    // next nearest-neighbor coupling
    if (sameMagnitude(dx, ref1))             coupling = -1.285354857;
    else if (sameMagnitude(dx, ref3))        coupling = -1.283593357;
    else if (sameMagnitude(dx, ref1+ref3))   coupling = -1.282563990;
  }  
  else if (sameMagnitude(dr, dr_ref3)) {    // next next nearest-neighbor coupling
    if (sameMagnitude(dx, ref2))             coupling = -1.017321802;
    else if (sameMagnitude(dx, ref3))        coupling = -1.019007206;
    else if (sameMagnitude(dx, ref2+ref3))   coupling = -1.017030791;
  }  
  else if (sameMagnitude(dr, dr_ref4)) {    // 4th nearest-neighbor coupling
    if (sameMagnitude(dx, dr_ref4))          coupling = 0.234128911;
    else if (sameMagnitude(dy, dr_ref4))     coupling = 0.235651112;
    else if (sameMagnitude(dz, dr_ref4))     coupling = 0.236648311;
  }
  else if (sameMagnitude(dr, dr_ref5)) {    // 5th nearest-neighbor coupling
    if (sameMagnitude(dx, ref3))               coupling = 0.343174662;
    else if (sameMagnitude(dx, ref1+ref3))     coupling = 0.339611863;
    else if (sameMagnitude(dx, ref2+ref3))     coupling = 0.342872395;
  }
  else if (sameMagnitude(dr, dr_ref6)) {    // 6th nearest-neighbor coupling
    if (sameMagnitude(dx, ref2))               coupling = 0.704679265;
    else if (sameMagnitude(dx, ref3))          coupling = 0.697818134;
    else if (sameMagnitude(dx, dr_ref4+ref1))  coupling = 0.703107429;
  }
  else if (sameMagnitude(dr, dr_ref7)) {    // 7th nearest-neighbor coupling
    if (sameMagnitude(dx, ref1))               coupling = 0.464462095;
    else if (sameMagnitude(dx, ref3))          coupling = 0.467291989;
    else if (sameMagnitude(dx, dr_ref4+ref2))  coupling = 0.468039834;
  }

  return -coupling;           // minus sign for antiferromagnetic coupling

}


/*
double CrystalStructure3D::assignExchangeCouplings_testing(double dx, double dy, double dz, double dr)
{

  const double dr_ref    {1.0};
  const double threshold {0.0001};
  double coupling        {0.0};

  // might replace it by the one in Utilities/CompareNumbers.hpp
  auto sameMagnitude = [=](double a, double b) -> bool {return fabs(fabs(a) - fabs(b)) < threshold; }; 

  if (sameMagnitude(dr, dr_ref)) {    // 4th nearest-neighbor coupling
    if (sameMagnitude(dx, dr_ref))          coupling = -1.0;
    else if (sameMagnitude(dy, dr_ref))     coupling = -1.0;
    else if (sameMagnitude(dz, dr_ref))     coupling = -1.0;
  }

  return coupling;

}
*/


// TODO: Ad hoc for our system for now.  All the reference dr's and coupling strengths should be read from input file (using readHamiltonianTerms()).
double CrystalStructure3D::assignDzyaloshinskiiMoriyaInteractions(double dz, double dr)
{
  
  const double dr_ref1   {0.613054};
  const double dr_ref2   {0.913759};
  const double dr_ref3   {0.957453};
  //const double dr_ref4   {1.0};
  //const double dr_ref5   {1.17296};
  //const double dr_ref6   {1.35461};
  //const double dr_ref7   {1.38445};

  const double ref1      {0.22956};
  const double ref2      {0.27044};
  const double ref3      {0.5};
  double coupling        {0.0};
  
  //auto sameMagnitude = [=](double a, double b) -> bool { return fabs(fabs(a) - fabs(b)) < threshold; };
  //auto sameSign      = [=](double a, double b) -> bool { return a * b > 0.0; };

  if (sameMagnitude(dr, dr_ref1)) {         // nearest-neighbor coupling
    if (sameMagnitude(dz, ref1))
      coupling = (dz > 0.0) ? -0.071692682 : 0.071692682;
    else if (sameMagnitude(dz, ref2))
      coupling = (dz > 0.0) ? -0.177878866 : 0.177878866;
    else if (sameMagnitude(dz, ref3))
      coupling = (dz > 0.0) ? 0.089626860 : -0.089626860;
  }
  else if (sameMagnitude(dr, dr_ref2)) {    // next nearest-neighbor coupling
    if (sameMagnitude(dz, ref1))
      coupling = (dz > 0.0) ? -0.047653359 : 0.047653359;
    else if (sameMagnitude(dz, ref3))
      coupling = (dz > 0.0) ? -0.069320477 : 0.069320477;
    else if (sameMagnitude(dz, ref1+ref3))
      coupling = (dz > 0.0) ? 0.055421967 : -0.055421967;
  }
  else if (sameMagnitude(dr, dr_ref3)) {    // next next nearest-neighbor coupling
    if (sameMagnitude(dz, ref2))
      coupling = (dz > 0.0) ? -0.021778776 : 0.021778776;
    else if (sameMagnitude(dz, ref3))
      coupling = (dz > 0.0) ? -0.099840235 : 0.099840235;
    else if (sameMagnitude(dz, ref2+ref3))
      coupling = (dz > 0.0) ? 0.024074747 : -0.024074747;
  }
/*
  else if (sameMagnitude(dr, dr_ref4)) {    // 4th nearest-neighbor coupling (sign to be confirmed)
    if (sameMagnitude(dx, dr_ref4))         // <-- sign problem!
      coupling = (dx > 0.0) ? -0.066329338 : 0.066329338;
    else if (sameMagnitude(dy, dr_ref4))    // <-- sign problem!
      coupling = (dy > 0.0) ? -0.046879378 : 0.046879378;
    else if (sameMagnitude(dz, dr_ref4))
      coupling = (dz > 0.0) ? 0.000973312 : -0.000973312;
  }
*/

  return -coupling;               // minus sign for antiferromagnetic coupling
  
}


ObservableType CrystalStructure3D::getExchangeInterations()
{

  ObservableType energy {0.0};

  for (unsigned int atomID=0; atomID<systemSize; atomID++) {
    for (auto neighbor : neighborList[atomID]) {
      energy += neighbor.J_ij * (spin[atomID].x * spin[neighbor.atomID].x + 
                                 spin[atomID].y * spin[neighbor.atomID].y + 
                                 spin[atomID].z * spin[neighbor.atomID].z);
    }
  } 

  return 0.5 * energy;          // the factor of 0.5 is for correcting double counting

}


ObservableType CrystalStructure3D::getDzyaloshinskiiMoriyaInterations()
{

  // Cross product between two spins
  // (spin[atomID].y * spin[neighbor.atomID].z - spin[atomID].z * spin[neighbor.atomID].y) +   // x-direction
  // (spin[atomID].z * spin[neighbor.atomID].x - spin[atomID].x * spin[neighbor.atomID].z) +   // y-direction
  // (spin[atomID].x * spin[neighbor.atomID].y - spin[atomID].y * spin[neighbor.atomID].x);    // z-direction

  ObservableType energy {0.0};

  for (unsigned int atomID=0; atomID<systemSize; atomID++) {
    for (auto neighbor : neighborList[atomID])
      energy += neighbor.D_ij * (spin[atomID].x * spin[neighbor.atomID].y - spin[atomID].y * spin[neighbor.atomID].x);    // z-direction only 
  }

  return 0.5 * energy;           // the factor of 0.5 is for correcting double counting

}


std::tuple<ObservableType, ObservableType, ObservableType, ObservableType> CrystalStructure3D::getMagnetization()
{

  ObservableType m1 {0.0};
  ObservableType m2 {0.0};
  ObservableType m3 {0.0};
  ObservableType m4 {0.0};

  for (unsigned int atomID=0; atomID<systemSize; atomID++) {
    m1 += spin[atomID].x;
    m2 += spin[atomID].y;
    m3 += spin[atomID].z;
  }

  m4 = sqrt(m1 * m1 + m2 * m2 + m3 * m3);

  return {m1, m2, m3, m4};

}


ObservableType CrystalStructure3D::getDifferenceInExchangeInterations()
{

  ObservableType energyChange {0.0};

  for (auto neighbor : neighborList[currentPosition]) {
    energyChange += neighbor.J_ij * ((spin[currentPosition].x - oldSpin.x) * spin[neighbor.atomID].x + 
                                     (spin[currentPosition].y - oldSpin.y) * spin[neighbor.atomID].y + 
                                     (spin[currentPosition].z - oldSpin.z) * spin[neighbor.atomID].z );
  }

  return energyChange;

}


ObservableType CrystalStructure3D::getDifferenceInDzyaloshinskiiMoriyaInterations()
{

  ObservableType energyChange {0.0};

  // z-direction only
  for (auto neighbor : neighborList[currentPosition]) {
    energyChange -= neighbor.D_ij * (oldSpin.x * spin[neighbor.atomID].y - oldSpin.y * spin[neighbor.atomID].x);
    energyChange += neighbor.D_ij * (spin[currentPosition].x * spin[neighbor.atomID].y - spin[currentPosition].y * spin[neighbor.atomID].x);       
  }

  return energyChange;

}
