#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the potential library
#
POTLIB = ../potlib_f90.a
MODDIR = ../../Modules
.SUFFIXES: .f90 .mod $(SUFFIXES)
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c -O4 
REL = 

SRCS = \
	v_e_laser.f90 

.f90.o:
	$(FC) $(FFLAGS) -I$(MODDIR) $<

all: $(POTLIB)
$(POTLIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(POTLIB) $(?:.f90=.o)

sources: $(SRCS)
	co -l RCS/*.f90
clean:
	rm *.o
