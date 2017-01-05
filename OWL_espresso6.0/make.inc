# make.inc.  Generated from make.inc.in by configure.

# compilation rules

.SUFFIXES :
.SUFFIXES : .o .c .f .f90

# most fortran compilers can directly preprocess c-like directives: use
# 	$(MPIF90) $(F90FLAGS) -c $<
# if explicit preprocessing by the C preprocessor is needed, use:
# 	$(CPP) $(CPPFLAGS) $< -o $*.F90
#	$(MPIF90) $(F90FLAGS) -c $*.F90 -o $*.o
# remember the tabulator in the first column !!!

.f90.o:
	$(MPIF90) $(F90FLAGS) -c $<

.cpp.o:
	$(CXX) $(CXXFLAGS)  -c $<

# .f.o and .c.o: do not modify

.f.o:
	$(F77) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS)  -c $<



# Top QE directory, not used in QE but useful for linking QE libs with plugins
# The following syntax should always point to TOPDIR:
#   $(dir $(abspath $(filter %make.sys,$(MAKEFILE_LIST))))

TOPDIR = /global/homes/s/syuk/3_Programs/6_OWL_6_0_Cori

# DFLAGS  = precompilation options (possible arguments to -D and -U)
#           used by the C compiler and preprocessor
# FDFLAGS = as DFLAGS, for the f90 compiler
# See include/defs.h.README for a list of options and their meaning
# With the exception of IBM xlf, FDFLAGS = $(DFLAGS)
# For IBM xlf, FDFLAGS is the same as DFLAGS with separating commas

# MANUAL_DFLAGS  = additional precompilation option(s), if desired
#                  BEWARE: it does not work for IBM xlf! Manually edit FDFLAGS
MANUAL_DFLAGS  =
#DFLAGS         =  -D__INTEL -D__DFTI -D__MPI -D__PARA -D__ELPA  -D__OPENMP -D__NON_BLOCKING_SCATTER
DFLAGS         =  -D__OPENMP -D__INTEL -D__FFTW -D__MPI -D__PARA -D__SCALAPACK
FDFLAGS        = $(DFLAGS) $(MANUAL_DFLAGS)

# IFLAGS = how to locate directories with *.h or *.f90 file to be included
#          typically -I../include -I/some/other/directory/
#          the latter contains .e.g. files needed by FFT libraries

IFLAGS         = -I../include

# MOD_FLAGS = flag used by f90 compiler to locate modules
# Each Makefile defines the list of needed modules in MODFLAGS

MOD_FLAG      = -I

# Compilers: fortran-90, fortran-77, C
# If a parallel compilation is desired, MPIF90 should be a fortran-90
# compiler that produces executables for parallel execution using MPI
# (such as for instance mpif90, mpf90, mpxlf90,...);
# otherwise, an ordinary fortran-90 compiler (f90, g95, xlf90, ifort,...)
# If you have a parallel machine but no suitable candidate for MPIF90,
# try to specify the directory containing "mpif.h" in IFLAGS
# and to specify the location of MPI libraries in MPI_LIBS

#MPIF90         = mpiifort
MPIF90		= ftn
#MPIF90		= tau_f90.sh

#F90		= ifort
F90           = ftn
#F90           = tau_f90.sh

#CC             = mpiicc
CC		= cc
CXX		= CC
#CC		= tau_cc.sh

#F77            = mpiifort
F77		= ftn
#F77		= tau_f90.sh

# C preprocessor and preprocessing flags - for explicit preprocessing,
# if needed (see the compilation rules above)
# preprocessing flags must include DFLAGS and IFLAGS

CPP            = cpp
CPPFLAGS       = -P -C -traditional $(DFLAGS) $(IFLAGS)

# compiler flags: C, F90, F77
# C flags must include DFLAGS and IFLAGS
# F90 flags must include MODFLAGS, IFLAGS, and FDFLAGS with appropriate syntax

CFLAGS         = -qopenmp -xCORE-AVX2 -g -O3 $(DFLAGS) $(IFLAGS)
CXXFLAGS       = -O3 -std=c++11 $(DFLAGS) $(IFLAGS)
F90FLAGS       = $(FFLAGS) -nomodule -qopenmp -fpp $(FDFLAGS) $(IFLAGS) $(MODFLAGS)
FFLAGS         = -debug -align array64byte -xCORE-AVX2 -O2 -assume byterecl -g -qopenmp -qopt-report=5 -qopt-report-phase=vec,openmp,loop

# compiler flags without optimization for fortran-77
# the latter is NEEDED to properly compile dlamch.f, used by lapack

FFLAGS_NOOPT   = -O0 -assume byterecl -g -traceback

# compiler flag needed by some compilers when the main program is not fortran
# Currently used for Yambo

FFLAGS_NOMAIN   = -nofor_main

# Linker, linker-specific flags (if any)
# Typically LD coincides with F90 or MPIF90, LD_LIBS is empty

#LD             = mpiifort
LD		= ftn
#LD             = tau_f90.sh
LDFLAGS        = -static-intel -qopenmp -nofor_main -cxxlib
LD_LIBS        = #-L$(MEMKIND_DIR)/lib -lmemkind

# External Libraries (if any) : blas, lapack, fft, MPI

# If you have nothing better, use the local copy :
# BLAS_LIBS = /your/path/to/espresso/BLAS/blas.a
# BLAS_LIBS_SWITCH = internal

BLAS_LIBS      =  
BLAS_LIBS_SWITCH = external

# If you have nothing better, use the local copy :
# LAPACK_LIBS = /your/path/to/espresso/lapack-3.2/lapack.a
# LAPACK_LIBS_SWITCH = internal
# For IBM machines with essl (-D__ESSL): load essl BEFORE lapack !
# remember that LAPACK_LIBS precedes BLAS_LIBS in loading order

LAPACK_LIBS    =   
LAPACK_LIBS_SWITCH = external

ELPA_LIBS_SWITCH = disabled
SCALAPACK_LIBS = 

# nothing needed here if the the internal copy of FFTW is compiled
# (needs -D__FFTW in DFLAGS)

FFT_LIBS       = 

# For parallel execution, the correct path to MPI libraries must
# be specified in MPI_LIBS (except for IBM if you use mpxlf)

MPI_LIBS       = 

# IBM-specific: MASS libraries, if available and if -D__MASS is defined in FDFLAGS

MASS_LIBS      = 

# ar command and flags - for most architectures: AR = ar, ARFLAGS = ruv

AR             = ar
ARFLAGS        = ruv

# ranlib command. If ranlib is not needed (it isn't in most cases) use
# RANLIB = echo

RANLIB         = ranlib

# all internal and external libraries - do not modify

FLIB_TARGETS   = all

LIBOBJS        = $(TOPDIR)/clib/clib.a $(TOPDIR)/iotk/src/libiotk.a
LIBS           = $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(FFT_LIBS) $(BLAS_LIBS) $(MPI_LIBS) $(MASS_LIBS) $(LD_LIBS)

# wget or curl - useful to download from network
WGET = wget -O

# Install directory - not currently used
PREFIX = /usr/common/software/espresso/5.4.0/intel
