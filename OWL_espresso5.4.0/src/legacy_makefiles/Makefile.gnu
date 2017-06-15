# Makefile for OWLPW

include ../../make.sys

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src \
	  $(MOD_FLAG)../../FFTXlib \
          $(MOD_FLAG)../../LAXlib \
	  $(MOD_FLAG)../../Modules \
	  $(MOD_FLAG)../../PW/src $(MOD_FLAG).
IFLAGS=

#location of needed libraries
LIBOBJS= ../../iotk/src/libiotk.a ../../clib/clib.a

COMMLIBS = \
Communications.o

OWLPWOBJS = \
WL_DFT_Interface.o

OWLPWLIBS = \
wl_getInfo_pwscf.o \
wl_do_pwscf.o \
wl_run_QE.o

MCOBJS = \
Histogram.o \
MCAlgorithms.o \
MCMoves.o

IOOBJS = \
InputOutput.o

QEMODS=../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a
PWOBJS= ../../PW/src/libpw.a 

TLDEPS=bindir mods libs liblapack libblas pw

all : tldeps owlpw.x

owlpw.x : libcomm.a $(OWLPWOBJS) $(MCOBJS) libowlpw.a $(IOOBJS) $(LIBOBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	$(OWLPWOBJS) $(MCOBJS) libowlpw.a libcomm.a $(IOOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin; ln -fs ../OWLPW/src/$@ . )

libcomm.a : $(COMMLIBS) 
	$(AR) $(ARFLAGS) $@ $?

libowlpw.a : $(OWLPWLIBS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/owlpw.x

veryclean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L *.cpp *.hpp
	- /bin/rm -f ../../bin/owlpw.x

include make.depend
