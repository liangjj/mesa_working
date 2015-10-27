#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the davidson eigenvalue and linear system routines
#
DVDLIB = dvdlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = f77
FFLAGS = -c 
GET = sccs get
REL = 

SRCS = \
	dvddat.f \
	frmres.f \
	honv.f \
	newvec.f \
	rdiag.f \
	reinit.f \
	rsdvd.f \
	tstovl.f 


.f.o:
	$(FC) $(FFLAGS) $<

all: $(DVDLIB)
$(DVDLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(DVDLIB) $(?:.f=.o)



