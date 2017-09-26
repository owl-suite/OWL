#ifndef TYPE_TRAITS_HPP
#define TYPE_TRAITS_HPP

#include <cmath>
#include <mpi.h>
#include <hdf5.h>


template<typename T>
class TypeTraits;


template<>
class TypeTraits<double>
{

public:
  inline static MPI_Datatype MPItype(void) { return MPI_DOUBLE; }
  inline static hid_t HDF5Type(void)       { return H5T_NATIVE_DOUBLE; }

};


template<>
class TypeTraits<int>
{

public:
  inline static MPI_Datatype MPItype(void) { return MPI_INT; }
  inline static hid_t HDF5Type(void)       { return H5T_NATIVE_INT; }

};


template<>
class TypeTraits<char>
{

public:
  inline static MPI_Datatype MPItype(void) { return MPI_CHAR; }
  inline static hid_t HDF5Type(void)       { return H5T_NATIVE_CHAR; }
};

#endif

