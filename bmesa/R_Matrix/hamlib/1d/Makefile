#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr hamiltonian routines
#
HAMLIB = ../hamlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c 
GET = sccs get
REL = 

SRCS = \
	drvh1.f \
	h1pac.f \
	prpac1.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(HAMLIB)
$(HAMLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(HAMLIB) $(?:.f=.o)



