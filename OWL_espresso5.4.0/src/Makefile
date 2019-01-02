### Define the current directory as the source directory
export SRC_DIR = $(shell pwd)

### Define variables
export BASIC_OWL_OBJS =                \
$(SRC_DIR)/Communications.o            \
$(SRC_DIR)/Heisenberg2D.o              \
$(SRC_DIR)/Ising2D.o                   \
$(SRC_DIR)/Histogram.o                 \
$(SRC_DIR)/RandomNumberGenerator.o     \
$(SRC_DIR)/MCAlgorithms.o              \
$(SRC_DIR)/WangLandauSampling.o        \
$(SRC_DIR)/ReplicaExchangeWangLandau.o \
$(SRC_DIR)/MulticanonicalSampling.o    \
$(SRC_DIR)/HistogramFreeMUCA.o         \
$(SRC_DIR)/InputOutput.o               \
$(SRC_DIR)/Globals.o                   \
$(SRC_DIR)/OWLMain.o

################################################

.PHONY: default all owl-qe clean

default : owl

all : owl owl-qe

owl : $(BASIC_OWL_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -o $@ $(BASIC_OWL_OBJS)

# Add "-DDRIVER_MODE_QE" to the CXXFLAGS to compile the QE-related part in those $BASIC_OWL_OBJS code 
owl-qe : export CXXFLAGS += -DDRIVER_MODE_QE
owl-qe :
	@if test -d QuantumEspresso ; then           \
	    cd QuantumEspresso && $(MAKE) owl-qe;    \
	fi

clean:
	rm -rf *.o *.a *.dSYM
	cd QuantumEspresso && $(MAKE) clean
	@for exe in owl owl-qe owl-feram owl-lsms; \
	do \
	    if test -f $$exe ; then rm $$exe; fi \
	done

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<

