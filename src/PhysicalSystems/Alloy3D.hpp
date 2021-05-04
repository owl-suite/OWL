#ifndef ALLOY3D_HPP
#define ALLOY3D_HPP

#include <filesystem>
#include <vector>
#include "CrystalBase.hpp"
#include "PhysicalSystemBase.hpp"
#include "Main/Globals.hpp"



//struct NeighboringAtom {
//  unsigned int atomID;
//  double       distance {0.0};                           // Distance from a reference atom
//  Element      atomic_species;
//};


class Alloy3D : public PhysicalSystem {

public :

  Alloy3D(const char* inputFile, int initial = 0); 
  ~Alloy3D();

  void writeConfiguration(int = 0, const char* = NULL)  override;
  void getObservables()                                 override;
  void doMCMove()                                       override;
  void acceptMCMove()                                   override;
  void rejectMCMove()                                   override;

  void getAdditionalObservables()                       override;

  //void buildMPIConfigurationType()                      override;
  //void readHamiltonianTerms(const char* inputFile);

  Lattice lattice;

private :

  // Model specific information
  unsigned int              numberOfElements {0};
  std::vector<Element>      elementTypes;
  std::vector<double>       composition;
  Matrix<double>            interactions;
  std::vector<unsigned int> numberOfAtomsForEachElement;
  std::vector<double>       neighborInteractionStrengths;

  // Overall configuration
  std::vector<Element> atom;

  // Old configuration
  unsigned int  currentPosition1;
  unsigned int  currentPosition2;
  Element       oldAtom1;
  Element       oldAtom2;
  
  bool firstTimeGetMeasures;

  // A matrix that stores nearest-neighhbor pairs in a configuration
  // One matrix for each nearest neighbor distance smaller than the nearestNeighborCutoff
  std::vector< Matrix<ObservableType> > nearestNeighborPairTypes;
  ObservableType                        idealEntropy      {0.0};
  ObservableType                        mutualInformation {0.0};

  // Private functions

  // Initialization:
  void         readCompositionInfo(const std::filesystem::path& mainInputFile);
  void         readAtomConfigFile(const std::filesystem::path& atomConfigFile);
  void         initializeAtomConfiguration(int initial = 0); 
  //double       assignExchangeCouplings(double dx, double dy, double dz, double dr);

  unsigned int getElementIndex(Element elem);

  // Hamiltonian measurements:
  void            getObservablesFromScratch();
  ObservableType  getIdealEntropy();
  ObservableType  getExchangeInteractions();
  ObservableType  getDifferenceInExchangeInteractions();
  void            getNearestNeighborPairTypes();
  ObservableType  getMutualInformation(unsigned int k = 0);
  ObservableType  getInformationEntropy();
 
};

#endif
