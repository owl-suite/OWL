#####################################################################
# Makefile for Oak-Ridge Wang-Landau (OWL)
# 
# Copyright (C) 2015-2018 Ying Wai Li, Oak Ridge National Laboratory
#
#####################################################################

### Define C++ compiler and flags for different machines
include architecture.inc

### Define directories
export MASTER_DIR = $(shell pwd)
export INCLUDE_PATH = 
export LIBRARY_PATH = 

################################################

.PHONY : default all owl owl-qe clean

default :
	@echo ' '
	@echo 'Please choose an option to compile the Oak-Ridge Wang-Landau (OWL) code with:'
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

all : owl owl-qe

owl :
	@if test -d src ; then           \
	    cd src && $(MAKE) owl;       \
	 fi

# Add "-DDRIVER_MODE_QE" to the CXXFLAGS to compile the QE-related part in the code
owl-qe : export CXXFLAGS += -DDRIVER_MODE_QE
owl-qe :
	@if test -d src ; then                            \
	    cd src && $(MAKE) owl-qe DRIVER_MODE_QE=1;    \
	 fi

clean :
	@if test -d bin ; then           \
	    cd bin && $(MAKE) clean;     \
	 fi
	@if test -d src ; then           \
	    cd src && $(MAKE) clean;     \
	 fi
