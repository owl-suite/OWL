# include <mpi.h>
# include "Matrix.hpp"          // Matrix class header
# include "MCMoves.hpp"
# include "WL_DFT_Interface.hpp"
# include "InputOutput.hpp"
# include "Histogram.hpp"
# include "Communications.hpp"


// These two should be derived from MonteCarloAlgorithm
void YingWaisCheck(int, int&);
void WangLandauSampling(int, int&, int);



class MonteCarloAlgorithm {

public:

  // Constructor
  MonteCarloAlgorithm();

  // Destructor
  ~MonteCarloAlgorithm();

protected:

  unsigned long int totalMCsteps;
  unsigned long int acceptedMoves;
  unsigned long int rejectedMoves;


private:


};

