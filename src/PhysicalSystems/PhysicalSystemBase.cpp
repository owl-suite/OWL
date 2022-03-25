#include <cassert>
#include <limits>
#include "PhysicalSystemBase.hpp"

// Specialized for spin models for now
void PhysicalSystem::calculateThermodynamics(std::vector<ObservableType> averagedObservables, std::vector<ObservableType> averagedObservablesSquared, double temperature)
{
  ObservableType specificHeat           {0.0};
  ObservableType magneticSusceptibility {0.0};
  ObservableType BinderCumulant         {0.0};
  unsigned int index = std::numeric_limits<unsigned int>::max();;

  for (unsigned int i=0; i < numObservables; i++) {

    //std::cout << "Observable name: " << observableName[i] << "\n";

    if (observableName[i] == "Total energy, E") {
      specificHeat = (averagedObservablesSquared[i] - averagedObservables[i] * averagedObservables[i]) / 
                     (systemSize * temperature * temperature);
      printf("   Specific heat, Cv          : %12.5f     (per site) \n", specificHeat);
      continue;
    }
    else if (observableName[i] == "Total absolute magnetization, |M|") {
      index = i;
      magneticSusceptibility = (averagedObservablesSquared[i] - averagedObservables[i] * averagedObservables[i]) / 
                               (systemSize * temperature);
      printf("   Magnetic susceptibility, \u03C7 : %12.5f     (per site) \n", magneticSusceptibility);
      continue;
    }
    else if (observableName[i] == "4th order magnetization, M^4") {
      assert (index < numObservables);
      BinderCumulant = 1.0 - averagedObservables[i] / (3.0 * averagedObservablesSquared[index] * averagedObservablesSquared[index]);
      printf("   Binder Cumulant, U4        : %12.5f \n",                BinderCumulant);
      continue;
    }

  }

  printf("\n");


  // Write results into an output file
  FILE* thermoFile;
  thermoFile = fopen("thermodynamics.dat", "w");

  fprintf(thermoFile, "# Thermodynamic quantities \n\n");
  fprintf(thermoFile, "Specific heat, Cv          : %12.5f     (per site) \n",      specificHeat);
  fprintf(thermoFile, "Magnetic susceptibility, \u03C7 : %12.5f     (per site) \n", magneticSusceptibility);
  fprintf(thermoFile, "Binder Cumulant, U4        : %12.5f \n",                     BinderCumulant);

  fclose(thermoFile);

}
