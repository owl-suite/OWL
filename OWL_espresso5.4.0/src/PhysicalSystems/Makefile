# Add the source directory to the include path
#export INCLUDE_PATH += -I $(SRC_DIR)

PHYS_SYS_OBJS = Ising2D.o       \
                Heisenberg2D.o

default : all

all : libPhysicalSystems.a

clean :
	rm -rf *.o *.a *.dSYM

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<

libPhysicalSystems.a : $(PHYS_SYS_OBJS)
	$(AR) -r $@ $(PHYS_SYS_OBJS)
