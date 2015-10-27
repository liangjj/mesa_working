#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the davidson eigenvalue and linear system routines
#
LINSYSLIB = ../linsyslib.a
AR = ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c -O5
GET = sccs get
REL = 

SRCS = \
	bsmall.f \
	lares.f \
	lavec.f \
	lindvd.f \
	lsolve.f \
	tstovl.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(LINSYSLIB)
$(LINSYSLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(LINSYSLIB) $(?:.f=.o)



