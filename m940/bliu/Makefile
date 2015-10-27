#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the davidson eigenvalue and linear system routines
#
BLIULIB = bliulib.a
GET = sccs get
REL = 

SRCS = \
	dvddat.f \
	frmres.f \
	honv.f \
	hsmall.f \
	newvec.f \
	rdiag.f \
	rdciv.f \
	rsdvd.f \
	tstovl.f 

all: $(BLIULIB)
$(BLIULIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(BLIULIB) $(?:.f=.o)

clean:
	rm -f *.o *.a