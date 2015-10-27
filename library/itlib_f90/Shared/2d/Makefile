#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr iterative routines
#
ITDVRLIB = ../../itdvrlib_f90.a
MODDIR = ../../../Modules
.SUFFIXES: .f90 .mod $(SUFFIXES)
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c -O5
GET = sccs get
REL = 

SRCS = \
	h22dvr.f90 \
	h2e.f90 \
	h22h0.f90 

.f90.o:
	$(FC) $(FFLAGS) -I$(MODDIR) $<

all: $(ITDVRLIB)
$(ITDVRLIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(ITDVRLIB) $(?:.f90=.o)
