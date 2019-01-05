# Makefile for OWL

# Add the source directory to the include path
#export INCLUDE_PATH += -I $(SRC_DIR)

OWL_MAIN_OBJS = Communications.o \
                InputOutput.o    \
                Globals.o        \
                OWLMain.o

default : owl

all : owl owl-qe

owl : $(OWL_MAIN_OBJS)
	#$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -o $@ $(OWL_MAIN_OBJS) $(BASIC_OWL_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -o $@ $(OWL_MAIN_OBJS) $(OWL_BASIC_LIBS)
	cp $@ $(MASTER_DIR)/bin

clean :
	rm -rf *.o *.a *.dSYM
	@for exe in owl owl-qe owl-feram owl-lsms; \
	do \
	    if test -f $$exe ; then rm $$exe; fi \
	done

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<
