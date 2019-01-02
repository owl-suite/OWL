export BASIC_OWL_OBJS = \
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


default : owl

all : owl owl-qe

owl : $(BASIC_OWL_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -o $@ $(BASIC_OWL_OBJS)

# TODO: Put "-DDRIVER_MODE_QE" in the CXXFLAGS to compile the QE-related part in those $BASIC_OWL_OBJS code after moving into QE subdirectory (Will the following work?)
#owl-qe : CXXFLAGS += -DDRIVER_MODE_QE

owl-qe : $(BASIC_OWL_OBJS) $(LIBOBJS) $(PWOBJS) $(QEMODS)
	$(MAKE) -f Makefile_QE owl-qe 
### TODO: To be moved into the QE subdirectory ###############################
#	if test -d qe_tmp ; then     \
#	    cd qe_tmp && $(MAKE) ;   \
#	fi
##############################################################################

clean:
	rm -rf *.o *.a *.dSYM
	@for exe in owl owl-qe owl-feram owl-lsms; \
	do \
	    if test -f $$exe ; then rm $$exe; fi \
	done

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<

