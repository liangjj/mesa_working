#   @(#)Makefile	5.1 11/6/94 
#
#    Makefile for the double precision integral library of MesaLib
#    
#    There are four routines partitioned off into SRCDP.
#    These need to be compiled as double-precision routines on both
#    32 and 64 bit machines.
#    
MESALIB = ../mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

DPFLAGS = -c
GET = sccs get
REL = 

SRCDP= \
	../dpintgrl/dfunc.f \
	../dpintgrl/dnode.f \
	../dpintgrl/droot.f \
	../dpintgrl/dsmit.f

.f.o:
	$(FC) $(DPFLAGS) $<

all: $(MESALIB)
$(MESALIB): $(SRCDP:.f=.o) 
	$(AR) $(ARFLAGS) $(MESALIB) $(?:.f=.o)
	$(RANLIB) $(MESALIB)

sources: $(SRCDP) 
$(SRCDP):
	$(GET) $(RELEASE) $@

link: sources
	rm -f ../../source/dpintgrl.f
	cat $(SRCDP) > ../../source/dpintgrl.f

print: link
	lpr ../../source/dpintgrl.f

clean:
	rm -f *.o
	rm -f core
