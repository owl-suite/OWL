#ifndef COMPARE_NUMBERS_HPP
#define COMPARE_NUMBERS_HPP

#include <cmath>

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


#endif