### Define the current directory as the source directory
export SRC_DIR = $(shell pwd)
export INCLUDE_PATH += -I $(SRC_DIR)

### Define variables
#export BASIC_OWL_OBJS =                \
#$(SRC_DIR)/Heisenberg2D.o              \
#$(SRC_DIR)/Ising2D.o                   \
#$(SRC_DIR)/RandomNumberGenerator.o     \
#$(SRC_DIR)/MCAlgorithms.o              \
#$(SRC_DIR)/Histogram.o                 \
#$(SRC_DIR)/WangLandauSampling.o        \
#$(SRC_DIR)/ReplicaExchangeWangLandau.o \
#$(SRC_DIR)/MulticanonicalSampling.o    \
#$(SRC_DIR)/HistogramFreeMUCA.o
#$(SRC_DIR)/Communications.o            \
#$(SRC_DIR)/InputOutput.o               \
#$(SRC_DIR)/Globals.o                   \
#$(SRC_DIR)/OWLMain.o

### Define libraries
export OWL_BASIC_LIBS = -L$(SRC_DIR)/Utilities -lUtilities                       \
                        -L$(SRC_DIR)/PhysicalSystems -lPhysicalSystems           \
                        -L$(SRC_DIR)/MonteCarloAlgorithms -lMonteCarloAlgorithms

################################################

.PHONY: default all owl-basics owl owl-qe clean

default : owl

all : owl owl-qe

owl-basics:
	cd Utilities && $(MAKE)
	cd PhysicalSystems && $(MAKE)
	cd MonteCarloAlgorithms && $(MAKE)

#owl : $(BASIC_OWL_OBJS)
#	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -o $@ $(BASIC_OWL_OBJS)
owl : owl-basics
	cd Main && $(MAKE) owl


# Add "-DDRIVER_MODE_QE" to the CXXFLAGS to compile the QE-related part in those $BASIC_OWL_OBJS code 
owl-qe : export CXXFLAGS += -DDRIVER_MODE_QE
owl-qe :
	@if test -d QuantumEspresso ; then           \
	    cd QuantumEspresso && $(MAKE) owl-qe;    \
	fi
# TODO:
#	cd Main && $(MAKE) owl-qe

clean:
	cd Main && $(MAKE) clean
	cd QuantumEspresso && $(MAKE) clean
	rm -rf *.o *.a *.dSYM
	@for exe in owl owl-qe owl-feram owl-lsms; \
	do \
	    if test -f $$exe ; then rm $$exe; fi \
	done

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<

