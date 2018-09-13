QE_OBJS =                   \
QuantumEspressoSystem.o     \
QEMCMoves.o                 \
owl_getInfo_pwscf.o         \
owl_do_pwscf.o              \
owl_run_QE.o
 

all : $(QE_OBJS)


clean :
	rm -f *.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<

%.o: %.f
	$(F77) $(FFLAGS) -c $< -o $@

