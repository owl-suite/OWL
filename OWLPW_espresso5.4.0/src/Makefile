# Makefile for OWLPW
QUANTUM_ESPRESSO_ROOT=./qe
include $(QUANTUM_ESPRESSO_ROOT)/make.sys

EXECUTABLE = owlpw

### C++ Compiler:
CXX = CC
### Compiler Flags:
#CXXFLAGS = -Wall -std=c++11 -g
CXXFLAGS = -O3 -Wall -std=c++11 -g

### Linking flags are compiler dependent
### For GNU compiler:
LINKFLAGS = $(FFLAGS_NOMAIN) -lstdc++

### For Intel compiler:
#LINKFLAGS = $(FFLAGS_NOMAIN) -cxxlib

.SUFFIXES :
.SUFFIXES : .o .c .f .f90 .cpp

.cpp.o:
	$(CXX) $(CXXFLAGS)  -c $<

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/iotk/src \
	  $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/FFTXlib  \
          $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/LAXlib   \
	  $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/Modules  \
	  $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/PW/src $(MOD_FLAG).
IFLAGS=

#location of needed libraries
LIBOBJS= $(QUANTUM_ESPRESSO_ROOT)/iotk/src/libiotk.a $(QUANTUM_ESPRESSO_ROOT)/clib/clib.a

COMMLIBS = \
Communications.o

OWLPWOBJS = \
QuantumEspressoSystem.o   \
Heisenberg2D.o            \
OWLMain.o

OWLPWLIBS = \
wl_getInfo_pwscf.o \
wl_do_pwscf.o \
wl_run_QE.o

MCOBJS = \
Histogram.o \
RandomNumberGenerator.o \
MCAlgorithms.o \
MCMoves.o

IOOBJS = \
InputOutput.o

QEMODS=$(QUANTUM_ESPRESSO_ROOT)/Modules/libqemod.a $(QUANTUM_ESPRESSO_ROOT)/FFTXlib/libqefft.a $(QUANTUM_ESPRESSO_ROOT)/LAXlib/libqela.a
PWOBJS= $(QUANTUM_ESPRESSO_ROOT)/PW/src/libpw.a

TLDEPS=bindir mods libs liblapack libblas pw

all : tldeps $(EXECUTABLE)

$(EXECUTABLE) : libcomm.a $(OWLPWOBJS) $(MCOBJS) libowlpw.a $(IOOBJS) $(LIBOBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) $(LINKFLAGS) -o $@ \
	$(OWLPWOBJS) $(MCOBJS) libowlpw.a libcomm.a $(IOOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd $(QUANTUM_ESPRESSO_ROOT)/bin; ln -fs $(CURDIR)/$@ . )

libcomm.a : $(COMMLIBS) 
	$(AR) $(ARFLAGS) $@ $?

libowlpw.a : $(OWLPWLIBS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd $(QUANTUM_ESPRESSO_ROOT) ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L $(EXECUTABLE)
	- /bin/rm -f $(QUANTUM_ESPRESSO_ROOT)/bin/$(EXECUTABLE)

veryclean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L *.cpp *.hpp $(EXECUTABLE)
	- /bin/rm -f $(QUANTUM_ESPRESSO_ROOT)/bin/$(EXECUTABLE)

include make.depend
