# Makefile for WLPW

include ../../make.sys

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src \
	  $(MOD_FLAG)../../FFTXlib \
	  $(MOD_FLAG)../../Modules \
	  $(MOD_FLAG)../../PW/src $(MOD_FLAG).
IFLAGS=

#location of needed libraries
LIBOBJS= ../../iotk/src/libiotk.a ../../flib/flib.a \
         ../../clib/clib.a ../../flib/ptools.a

WLPWOBJS = \
WL_DFT_Interface.o

WLPWLIBS = \
wl_getInfo_pwscf.o \
wl_do_pwscf.o \
wl_run_QE.o

MCLIBS = \
MCAlgorithms.o \
MCMoves.o

QEMODS=../../Modules/libqemod.a ../../FFTXlib/libqefft.a
PWOBJS= ../../PW/src/libpw.a

TLDEPS=bindir mods libs liblapack libblas pw

all : tldeps wlpw.x

wlpw.x : $(WLPWOBJS) libwlpw.a libmc.a $(LIBOBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ -nofor_main -cxxlib \
	$(WLPWOBJS) libwlpw.a libmc.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin; ln -fs ../WLPW/src/$@ . )

libwlpw.a : $(WLPWLIBS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

libmc.a : $(MCLIBS)
	$(AR) $(ARFLAGS) $@ $?

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/wlpw.x

include make.depend
