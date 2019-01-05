# Add the source directory to the include path
#export INCLUDE_PATH += -I $(SRC_DIR)

UTIL_OBJS = RandomNumberGenerator.o

default : all

all : libUtilities.a

clean :
	rm -rf *.o *.a *.dSYM

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<

libUtilities.a : $(UTIL_OBJS)
	$(AR) -r $@ $(UTIL_OBJS)
