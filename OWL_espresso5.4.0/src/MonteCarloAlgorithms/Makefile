# Add the source directory to the include path
#export INCLUDE_PATH += -I $(SRC_DIR)

MONTE_CARLO_OBJS = MCAlgorithms.o               \
                   Histogram.o                  \
                   MulticanonicalSampling.o     \
                   WangLandauSampling.o         \
                   ReplicaExchangeWangLandau.o  \
                   HistogramFreeMUCA.o

default : all

all : libMonteCarloAlgorithms.a

clean :
	rm -rf *.o *.a *.dSYM

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<

libMonteCarloAlgorithms.a : $(MONTE_CARLO_OBJS)
	$(AR) -r $@ $(MONTE_CARLO_OBJS)
