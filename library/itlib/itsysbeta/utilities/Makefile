#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the davidson eigenvalue and linear system routines
#
ITERLIB = ../iterlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c -O5
GET = sccs get
REL = 

SRCS = \
	abschm.f \
	copy.f \
	gschmt.f \
	honv.f \
	itoc.f \
	prntfm.f \
	prntrm.f \
	rzero.f 


.f.o:
	$(FC) $(FFLAGS) $<

all: $(ITERLIB)
$(ITERLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(ITERLIB) $(?:.f=.o)




