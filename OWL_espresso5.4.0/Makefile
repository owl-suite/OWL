#####################################################################
# Makefile for Oak-ridge Wang-Landau (OWL)
# 
# Copyright (C) 2015-2018 Ying Wai Li, Oak Ridge National Laboratory
#
#####################################################################

### Define C++ compiler and flags for different machines
#include architecture/mac_osx
#include architecture/olcf_titan
include architecture/olcf_eos
#include architecture/nersc

### Directory of Include-files:
export INCLUDE_PATH = -I. 
export LIBRARY_PATH = 

################################################

default :
	@echo ' '
	@echo 'Please choose an option to compile the Oak-ridge Wang-Landau (OWL) code with:'
	@echo 'make [options]'
	@echo ' '
	@echo 'Available package options:'
	@echo ' [options]             [descriptions]'
	@echo '  owl                   OWL stand-alone mode'
	@echo '  owl-qe                OWL interfaced with Quantum Espresso'
	@echo '  owl-lsms              OWL interfaced with Locally Self-consistent Multiple Scattering (LSMS)'
	@echo '  owl-feram             OWL interfaced with FERAM'
	@echo '  all                   Compile all owl, owl-qe, owl-lsms, owl-feram'
	@echo ' '
	@echo 'Available operation options:'
	@echo '  clean                 Remove OWL executables and objects'
	@echo '  cleanall (tentative)  Remove executables and objects for OWL and external codes'
	@echo ' '

owl :
	@if test -d src ; then           \
	    cd src && $(MAKE) owl;       \
	 fi

owl-qe :
	@if test -d src ; then                            \
	    cd src && $(MAKE) owl-qe DRIVER_MODE_QE=1     \
	 fi

all : owl

clean :
	@if test -d src ; then           \
	    cd src && $(MAKE) clean;     \
	 fi

