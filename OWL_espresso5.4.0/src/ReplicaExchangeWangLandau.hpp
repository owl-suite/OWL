#ifndef REPLICA_EXCHANGE_WANG_LANDAU_HPP
#define REPLICA_EXCHANGE_WANG_LANDAU_HPP

#include "MCAlgorithms.hpp"
#include "Histogram.hpp"
#include "Communications.hpp"


typedef int WalkerIDType;
typedef int WindowIDType;


class ReplicaExchangeWangLandau : public MonteCarloAlgorithm {

public :
  
  ReplicaExchangeWangLandau(MPICommunicator PhySystemComm, MPICommunicator MCAlgorithmComm);
  ~ReplicaExchangeWangLandau();

  void run(PhysicalSystem* physical_system);

private :

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

  int swapDirection;

  int upExchanges;
  int downExchanges;


  // Private member functions:
  void assignSwapPartner();
  void swapEnergy(double &energyForSwap);
  bool determineAcceptance(double myDOSRatio);
  void swapConfiguration(double configForSwap[], int numElements);
 
  void readREWLInputFile(const char* fileName); 

};


#endif
