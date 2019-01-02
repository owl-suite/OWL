BASIC_OWL_OBJS =            \
Communications.o            \
Heisenberg2D.o              \
Ising2D.o                   \
Histogram.o                 \
RandomNumberGenerator.o     \
MCAlgorithms.o              \
WangLandauSampling.o        \
ReplicaExchangeWangLandau.o \
MulticanonicalSampling.o    \
HistogramFreeMUCA.o         \
InputOutput.o               \
Globals.o                   \
OWLMain.o

### To be moved into the QE subdirectory #####################################

QUANTUM_ESPRESSO_ROOT=./qe
include $(QUANTUM_ESPRESSO_ROOT)/make.sys

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/iotk/src \
          $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/FFTXlib  \
          $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/LAXlib   \
          $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/Modules  \
          $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/PW/src $(MOD_FLAG).
#IFLAGS=

#location of needed libraries
LIBOBJS= $(QUANTUM_ESPRESSO_ROOT)/iotk/src/libiotk.a $(QUANTUM_ESPRESSO_ROOT)/clib/clib.a

QEMODS=$(QUANTUM_ESPRESSO_ROOT)/Modules/libqemod.a $(QUANTUM_ESPRESSO_ROOT)/FFTXlib/libqefft.a $(QUANTUM_ESPRESSO_ROOT)/LAXlib/libqela.a
PWOBJS= $(QUANTUM_ESPRESSO_ROOT)/PW/src/libpw.a

#TLDEPS=bindir mods libs liblapack libblas pw



# QE_OBJS should be compiled inside the QE subdirectory
QE_OBJS =                   \
QuantumEspressoSystem.o     \
QEMCMoves.o

OWLLIBS = \
owl_getInfo_pwscf.o         \
owl_do_pwscf.o              \
owl_run_QE.o

##############################################################################

default : owl

all : owl owl-qe

#owl : OWLMain.cpp $(BASIC_OWL_OBJS)
owl : $(BASIC_OWL_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -o $@ $(BASIC_OWL_OBJS)

#owl-qe : OWLMain.cpp $(BASIC_OWL_OBJS) $(QE_OBJS) libowl.a $(LIBOBJS) $(PWOBJS) $(QEMODS)
owl-qe : $(BASIC_OWL_OBJS) $(QE_OBJS) libowl.a $(LIBOBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) $(LINKFLAGS) $(FFLAGS_NOMAIN) -o $@ $(BASIC_OWL_OBJS) $(QE_OBJS) libowl.a $(PWOBJS) $(QEMODS) $(LIBOBJS) -L/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin $(LIBS)
# -L/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin
#owl-qe : OWLMain.cpp $(BASIC_OWL_OBJS) $(QE_OBJS)
### To be moved into the QE subdirectory #####################################
#	if test -d qe_tmp ; then     \
#	    cd qe_tmp && $(MAKE) ;   \
#	fi
##############################################################################
#	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -DDRIVER_MODE_QE -o owl-qe OWLMain.cpp $(BASIC_OWL_OBJS) $(QE_OBJS)

clean:
	rm -rf *.o *.a *.dSYM
	@for exe in owl owl-qe owl-feram owl-lsms; \
	do \
	    if test -f $$exe ; then rm $$exe; fi \
	done

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<
#	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c $*.cpp 

### To be moved into the QE subdirectory #####################################
#%.o: %.f
#	$(F77) $(FFLAGS) -c $< -o $@
libowl.a : $(OWLLIBS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

include make.depend
##############################################################################
