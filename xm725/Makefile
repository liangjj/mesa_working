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
TCGMSGLIB = /usr4/russo/src/pnl_stuff/tcgmsg/libtcgmsg.$(OPSYS).a
MDLIB =
GET = sccs get
REL = 

SRCS = \
	../xm725/d1xc.f \
	../xm725/dxcint.f \
	../xm725/dxclos.f \
	../xm725/dxopen.f \
	../xm725/fmaccum.f \
	../xm725/lebedev.f \
	../xm725/m725.f \
	../xm725/mkatmg.f \
	../xm725/pbasis.f \
	../xm725/plnkerr.f \
	../xm725/rhomax.f \
	../xm725/sg1.f \
	../xm725/sphere.f \
	../xm725/pm725.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m725 $(BINDIR)/m725.slave

../xm725/getscm.small.o: ../xm725/getscm.small.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) $<

../xm725/getscm.big.o: ../xm725/getscm.big.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) $<


$(BINDIR)/m725: $(SRCS:.f=.o) $(MESALIB) ../xm725/getscm.big.o
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../xm725/getscm.big.o $(MESALIB) $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m725

$(BINDIR)/m725.slave: $(SRCS:.f=.o) $(MESALIB) ../xm725/getscm.small.o
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../xm725/getscm.small.o $(MESALIB) $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m725.slave

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


