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
InputOutput.o

QE_OBJS =                   \
QuantumEspressoSystem.o     \
QEMCMoves.o                 \
owl_getInfo_pwscf.o         \
owl_do_pwscf.o              \
owl_run_QE.o


default : owl

all : owl owl-qe

owl : OWLMain.cpp $(BASIC_OWL_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -o owl OWLMain.cpp $(BASIC_OWL_OBJS)

owl-qe : OWLMain.cpp $(BASIC_OWL_OBJS) $(QE_OBJS)
	if test -d qe_tmp ; then     \
	    cd qe_tmp && $(MAKE) ;   \
	fi
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -DDRIVER_MODE_QE -o owl-qe OWLMain.cpp $(BASIC_OWL_OBJS) $(QE_OBJS)

clean:
	rm -rf *.o *.dSYM
	@for exe in owl owl-qe owl-feram owl-lsms; \
	do \
	    if test -f $$exe ; then rm $$exe; fi \
	done

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<
#	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c $*.cpp 
