#ifndef COMPARE_NUMBERS_HPP
#define COMPARE_NUMBERS_HPP

#include <cmath>
#include <vector>
#include <algorithm>

const double threshold {0.000001};

template <typename numberType>
inline bool sameMagnitude(numberType a, numberType b) { return fabs(fabs(a) - fabs(b)) < threshold; };


template <typename numberType>
inline bool sameSign(numberType a, numberType b) { return a * b > 0.0; };


template <typename numberType>
bool isFoundInVector(numberType a, std::vector<numberType> vec) {

  bool isFound {false};

  for (auto i : vec)
    if (sameMagnitude(a, i) && sameSign(a, i) ) {
        isFound = true;
        break;
    }

  return isFound;

}


// Returns the index of the vector at which the element is the same as the number in comparison
template <typename numberType>
unsigned int getVectorIndex(numberType a, std::vector<numberType> vec) {

  unsigned int index {0};
  for (unsigned int i=0; i<vec.size(); i++)
    if (sameMagnitude(a, vec[i]) && sameSign(a, vec[i])) {
      index = i;
      break;
    }

  return index;

}


// Returns the index of a sorted vector at which the element is greater than the number in comparison
template <typename numberType>
unsigned int getVectorUpperIndex(numberType a, std::vector<numberType> vec) {

  unsigned int index {0};
  for (unsigned int i=0; i<vec.size(); i++)
    if (a < vec[i]) {
      index = i;
      break;
    }

  return index;

}

#endif