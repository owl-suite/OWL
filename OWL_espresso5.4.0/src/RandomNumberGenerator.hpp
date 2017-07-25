#ifndef RANDOM_NUMBER_GENERATOR_HPP
#define RANDOM_NUMBER_GENERATOR_HPP

#include <random>
#include "Communications.hpp"

// An RNG for global use
extern std::mt19937 rng;
//extern int RngSeed;

void initializeRandomNumberGenerator(MPICommunicator, int = -1);


// Returns a random number between [-0.5,0.5]
inline double getRandomNumber()
{
  return (double(rng()) / double(rng.max()) - 0.5);
}

// Returns a random number between [0,1]
inline double getRandomNumber2()
{
  return (double(rng()) / double(rng.max()));
}

#endif

