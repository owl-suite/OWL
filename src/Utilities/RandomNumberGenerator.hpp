#ifndef RANDOM_NUMBER_GENERATOR_HPP
#define RANDOM_NUMBER_GENERATOR_HPP

#include <random>
#include "Main/Communications.hpp"

// TODO: allow for different choice of random number generators

// An RNG for global use
extern std::mt19937 rng_engine;
//extern int RngSeed;

// Uniform distribution between (-0.5, 0.5)
extern std::uniform_real_distribution<double> distribution1; 

// Uniform distribution between (0,1)
extern std::uniform_real_distribution<double> distribution2; 

// Uniform integer distribution between (0,?)
extern std::uniform_int_distribution<int> distribution_int; 

void initializeRandomNumberGenerator(MPICommunicator, int = -1);


// Returns a random number between (-0.5, 0.5)
inline double getRandomNumber()
{
  return distribution1(rng_engine);
}

// Returns a random number between (0.0, 1.0)
inline double getRandomNumber2()
{
  return distribution2(rng_engine);
}

// Returns an integral random number between (0,?)
inline int getIntRandomNumber()
{
  return distribution_int(rng_engine);
}

#endif

