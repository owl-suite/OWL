### Makefile for Main/ in OWL

### Define OWL's main objects
OWL_MAIN_OBJS = Communications.o \
                InputOutput.o    \
                Globals.o        \
                OWLMain.o

default : owl

all : libMain.a owl owl-qe

owl : libMain.a
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -o $@ $(OWL_LIBS)
	cp $@ $(MASTER_DIR)/bin

owl-qe : libMain.a

clean :
	rm -rf *.o *.a *.dSYM
	if test -f owl ; then rm owl; fi

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<

libMain.a : $(OWL_MAIN_OBJS)
	$(AR) -r $@ $(OWL_MAIN_OBJS)
	@echo " "
	