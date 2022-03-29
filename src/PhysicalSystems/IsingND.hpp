#ifndef ISINGND_HPP
#define ISINGND_HPP

#include <filesystem>
#include "PhysicalSystemBase.hpp"

typedef int          IsingSpinDirection;
typedef unsigned int Coordinates;


class IsingND : public PhysicalSystem {

public :

  IsingND(const char* spinConfigFile = "config_initial.dat", int = 0); 
  ~IsingND();

  //void readCommandLineOptions()                         override;
  void writeConfiguration(int = 0, const char* = NULL)  override;
  void getObservables()                                 override;
  void doMCMove()                                       override;
  void acceptMCMove()                                   override;
  void rejectMCMove()                                   override;

  void buildMPIConfigurationType();

private :

  unsigned int Size;
  unsigned int dimension;

  // Old configuration
  std::vector<unsigned int>  currentPosition;
  IsingSpinDirection         oldSpin;
  indexType                  currentIndex;

  // A 1D vector to store the spin configuration (a C-style array for MPI to operate on)
  IsingSpinDirection* spin;

  // Offsets to facilitate conversion between index and coordinates
  std::vector<indexType>     offsets;

  // Initialization:
  void      readSpinConfigFile(const std::filesystem::path& spinConfigFile);
  void      initializeSpinConfiguration(int initial = 0);
  indexType getIndexFromCoordinates(unsigned int dim, ...);                           // ... is the ND coordinates of a spin(x1, x2, x3, x4 ...)
  indexType getIndexFromCoordinates(std::vector<unsigned int> coords);
  void      getCoordinatesFromIndex(indexType index, std::vector<unsigned int>& coords);
  void      calculateOffsets();

};

#endif
