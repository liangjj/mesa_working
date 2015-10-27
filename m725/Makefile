# %W%   %G%
#
# Makefile for M725
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
TCGMSGLIB = $(MESA_LIBS)/tcgmsg/libtcgmsg.$(OPSYS).a
MDLIB =
GET = sccs get
REL = 

SRCS = \
	../m725/d1xc.f \
	../m725/dxcint.f \
	../m725/dxclos.f \
	../m725/dxopen.f \
	../m725/fmaccum.f \
	../m725/lebedev.f \
	../m725/m725.f \
	../m725/mkatmg.f \
	../m725/pbasis.f \
	../m725/plnkerr.f \
	../m725/rhomax.f \
	../m725/sg1.f \
	../m725/sphere.f \
	../m725/pm725.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m725 $(BINDIR)/m725.slave

../m725/getscm.small.o: ../m725/getscm.small.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) ../m725/getscm.small.$(OPSYS).f

../m725/getscm.big.o: ../m725/getscm.big.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) ../m725/getscm.big.$(OPSYS).f


$(BINDIR)/m725: $(SRCS:.f=.o) $(MESALIB) ../m725/getscm.big.o
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../m725/getscm.big.o $(MESALIB) $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m725

$(BINDIR)/m725.slave: $(SRCS:.f=.o) $(MESALIB) ../m725/getscm.small.o
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../m725/getscm.small.o $(MESALIB) $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m725.slave

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m725.f
	cat $(SRCS) > ../source/m725.f

print: link
	lpr ../source/m725.f

clean:
	rm -f *.o core


