#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the davidson eigenvalue and linear system routines
#
BLIULIB = ../bliulib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c -O5
GET = sccs get
REL = 

SRCS = \
	dvddat.f \
	frmres.f \
	hsmall.f \
	newvec.f \
	rdiag.f \
	rdciv.f \
	rsdvd.f \
	tstovl.f 


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BLIULIB)
$(BLIULIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(BLIULIB) $(?:.f=.o)




