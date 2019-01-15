# Makefile to clean binaries

default : all

all : clean

clean : 
	@for exe in owl owl-qe owl-feram owl-lsms; \
	do \
	    if test -f $$exe ; then rm $$exe; fi \
	done

