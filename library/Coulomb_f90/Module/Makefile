#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
include ../Makefile.inc
.SUFFIXES: .f90 .mod $(SUFFIXES)

SRCS = \
	Coulomb_Module.f90

all: $(COUL_LIB) 
$(COUL_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(COUL_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod
moddir:
	mv *.mod ../.mod_dir/.
