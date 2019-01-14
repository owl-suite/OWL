### Define the current directory as the source directory
export SRC_DIR = $(shell pwd)
export INCLUDE_PATH += -I $(SRC_DIR)

### Define libraries
export OWL_LIBS = -L$(SRC_DIR)/Utilities -lUtilities                       \
                  -L$(SRC_DIR)/PhysicalSystems -lPhysicalSystems           \
                  -L$(SRC_DIR)/MonteCarloAlgorithms -lMonteCarloAlgorithms \
                  -L$(SRC_DIR)/Main -lMain

################################################

.PHONY: default all owl-basics owl owl-qe clean

default : owl

all : owl owl-qe

owl-basics:
	cd Utilities && $(MAKE)
	cd PhysicalSystems && $(MAKE)
	cd MonteCarloAlgorithms && $(MAKE)

owl : owl-basics
	cd Main && $(MAKE) owl

# Add "-DDRIVER_MODE_QE" to the CXXFLAGS to compile the QE-related part in those $BASIC_OWL_OBJS code 
# Might want to move it to the Makefile up a level (the master one)
#owl-qe : export CXXFLAGS += -DDRIVER_MODE_QE
owl-qe : owl-basics
	cd Main && $(MAKE) owl-qe
	cd QuantumEspresso && $(MAKE) owl-qe

clean:
	cd Main && $(MAKE) clean
	cd Utilities && $(MAKE) clean
	cd PhysicalSystems && $(MAKE) clean 
	cd MonteCarloAlgorithms && $(MAKE) clean
	cd QuantumEspresso && $(MAKE) clean

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<

