# Makefile for OWL-QE

# $(LD) $(LDFLAGS) $(IFLAGS) $(LIBS) $(MPIF90) $(F90FLAGS) are in QE's make.sys
# $(LINKFLAGS) is in OWL's architecture file

### Quantum Espresso related stuff
###############################################################

# TODO: Test if "qe" links to a valid QE directory
# TODO: Download and compile QE if it is not there yet
QUANTUM_ESPRESSO_ROOT=./qe
include $(QUANTUM_ESPRESSO_ROOT)/make.sys
include make.depend

# location of needed modules and included files (if any)
MODFLAGS = $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/iotk/src \
           $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/FFTXlib  \
           $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/LAXlib   \
           $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/Modules  \
           $(MOD_FLAG)$(QUANTUM_ESPRESSO_ROOT)/PW/src $(MOD_FLAG).
IFLAGS =

#location of needed libraries
LIBOBJS = $(QUANTUM_ESPRESSO_ROOT)/iotk/src/libiotk.a $(QUANTUM_ESPRESSO_ROOT)/clib/clib.a

QEMODS = $(QUANTUM_ESPRESSO_ROOT)/Modules/libqemod.a $(QUANTUM_ESPRESSO_ROOT)/FFTXlib/libqefft.a $(QUANTUM_ESPRESSO_ROOT)/LAXlib/libqela.a
PWOBJS = $(QUANTUM_ESPRESSO_ROOT)/PW/src/libpw.a

#TLDEPS=bindir mods libs liblapack libblas pw


###############################################################

# Add the source directory to the include path
#export INCLUDE_PATH += -I $(SRC_DIR)

QE_OBJS =                   \
QuantumEspressoSystem.o     \
QEMCMoves.o                 \
owl_getInfo_pwscf.o         \
owl_do_pwscf.o              \
owl_run_QE.o

default: owl-qe

all : owl-qe

# TODO: Remove the MKL library link
# -L/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin
owl-qe : $(BASIC_OWL_OBJS) $(LIBOBJS) $(PWOBJS) $(QEMODS) $(QE_OBJS) $(OWLLIBS)
	$(LD) $(LDFLAGS) $(LINKFLAGS) $(FFLAGS_NOMAIN) -o $@ $(BASIC_OWL_OBJS) $(QE_OBJS) $(OWLLIBS) $(PWOBJS) $(QEMODS) $(LIBOBJS) -L/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin $(LIBS)
	cp $@ $(MASTER_DIR)/bin

clean:
	rm -rf *.o *.a *.dSYM owl-qe

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<

%.o: %.f90
	$(MPIF90) $(F90FLAGS) -c $< -o $@

