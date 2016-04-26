#include "InputOutput.hpp"

void writeEnergyFile(char fileName[], double energy)
{ 
  FILE *energy_file;
  energy_file = fopen(fileName, "a");
  fprintf(energy_file, "%f", energy);
  fclose(energy_file); 
} 

