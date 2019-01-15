# Makefile for QuantumEspresso/ in OWL

# $(LD) $(LDFLAGS) $(LIBS) $(MPIF90) $(F90FLAGS) are in QE's make.sys
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

#location of needed libraries
LIBOBJS = $(QUANTUM_ESPRESSO_ROOT)/iotk/src/libiotk.a $(QUANTUM_ESPRESSO_ROOT)/clib/clib.a

QEMODS = $(QUANTUM_ESPRESSO_ROOT)/Modules/libqemod.a $(QUANTUM_ESPRESSO_ROOT)/FFTXlib/libqefft.a $(QUANTUM_ESPRESSO_ROOT)/LAXlib/libqela.a
PWOBJS = $(QUANTUM_ESPRESSO_ROOT)/PW/src/libpw.a

# TODO: Check why it was introduced here originally
# YingWai: this might be a good starting point to compile QE if not present
#TLDEPS=bindir mods libs liblapack libblas pw


###############################################################


### Define OWL's QuantumEspresso interface objects
OWL_QE_OBJS = QuantumEspressoSystem.o     \
              QEMCMoves.o                 \
              owl_getInfo_pwscf.o         \
              owl_do_pwscf.o              \
              owl_run_QE.o

### Add the QuantumEspresso library to link line
export OWL_LIBS += -L$(SRC_DIR)/PhysicalSystems/QuantumEspresso -lQuantumEspresso

#owl-qe : tldeps libQuantumEspresso.a
owl-qe : libQuantumEspresso.a
	$(LD) $(LDFLAGS) $(LINKFLAGS) $(FFLAGS_NOMAIN) -o $@ $(OWL_LIBS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	cp $@ $(MASTER_DIR)/bin

clean:
	rm -rf *.o *.a *.dSYM 
	if test -f owl-qe ; then rm owl-qe; fi

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<

%.o: %.f90
	$(MPIF90) $(F90FLAGS) -c $< -o $@

libQuantumEspresso.a : $(OWL_QE_OBJS) 
	$(AR) -r $@ $(OWL_QE_OBJS)
	@echo " "

#tldeps :
#	if test -n "$(TLDEPS)" ; then \
#	( cd $(QUANTUM_ESPRESSO_ROOT) ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi
