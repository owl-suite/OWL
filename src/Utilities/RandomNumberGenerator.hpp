#ifndef RANDOM_NUMBER_GENERATOR_HPP
#define RANDOM_NUMBER_GENERATOR_HPP

#include <random>
#include "Main/Communications.hpp"

// TODO: allow for different choice of random number generators

// An RNG for global use
extern std::mt19937 rng;
//extern int RngSeed;

// Uniform distribution between (-0.5, 0.5)
extern std::uniform_real_distribution<double> distribution1; 

// Uniform distribution between (0,1)
extern std::uniform_real_distribution<double> distribution2; 


void initializeRandomNumberGenerator(MPICommunicator, int = -1);


// Returns a random number between (-0.5, 0.5)
inline double getRandomNumber()
{
  return distribution1(rng);
}

// Returns a random number between (0, 1)
inline double getRandomNumber2()
{
  return distribution2(rng);
}

#endif

