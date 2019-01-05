#ifndef REPLICA_EXCHANGE_WANG_LANDAU_HPP
#define REPLICA_EXCHANGE_WANG_LANDAU_HPP

#include "MCAlgorithms.hpp"
#include "Histogram.hpp"
#include "Main/Communications.hpp"

/*
  ReplicaExchange class:

  This class implements replica-exchange Wang-Landau sampling (REWL). 
  Reference: Phys. Rev. Lett. 110, 210603 (2013)
*/


// TO DO: distinguish WalkerType and WindowType
typedef int WalkerIDType;
typedef int WindowIDype;


class ReplicaExchangeWangLandau : public MonteCarloAlgorithm {

public :
  
  ReplicaExchangeWangLandau(PhysicalSystem* ps, MPICommunicator PhySystemComm, MPICommunicator MCAlgorithmComm);
  ~ReplicaExchangeWangLandau();

  void run();

private :

  PhysicalSystem* physical_system;
  Histogram h;

  /// YingWai's note: (Sep 17, 2017)
  /// Is it better to have the random number generator here?
  //std::mt19937 rng;
  //std::uniform_real_distribution<double> rnd;

  MPICommunicator PhysicalSystemComm;
  MPICommunicator REWLComm;
  int numWalkers;
  int numWindows;
  int numWalkersPerWindow;
  double overlap;

  int walkerID;
  int myWindow;
  int partnerID;
  int partnerWindow;

  int replicaExchangeInterval;
  int swapDirection;

  int upExchanges;
  int downExchanges;

  double MaxModFactor;


  // Private member functions:
  bool replicaExchange();                                    // Behave like a doMCMove; only propose a new energy and a new configuration
  void assignSwapPartner();
  void exchangeEnergy(ObservableType &energyForSwap);        // Swap energy with partner; energyForSwap will be overwritten
  bool determineAcceptance(double localDOSRatio);
  
  //template <typename T>
  //void exchangeConfiguration(T configForSwap[], int numElements);  
  void exchangeConfiguration(void* ptrToConfig, int numElements, MPI_Datatype MPI_config_type);

  void getMaxModFactor();

  void readREWLInputFile(const char* fileName); 

};


/*
template <typename T>
void ReplicaExchangeWangLandau::exchangeConfiguration(T configForSwap[], int numElements)
{

  if (partnerID != -1)
    REWLComm.swapVector(configForSwap, numElements, partnerID);

}
*/


#endif
