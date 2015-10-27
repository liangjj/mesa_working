#  @(#)Makefile	5.1 11/6/94
#
# Makefile for m613
#
FC = f77
FFLAGS = -c 

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB = 
GET = sccs get
REL = 

SRCS = \
	../m613/bakwrd.f \
	../m613/cotes.f \
	../m613/forwrd.f \
	../m613/ftogrid.f \
	../m613/ftoylm.f \
	../m613/gofish.f \
	../m613/interv.f \
	../m613/laplace.f \
	../m613/m613.f \
	../m613/ncwts.f \
	../m613/necote.f \
	../m613/newton.f \
	../m613/pm613.f \
	../m613/rsolver.f \
	../m613/rules.f \
	../m613/spline3.f \
	../m613/splint.f \
	../m613/tstylm.f \
	../m613/v0.f \
	../m613/v0atr.f \
	../m613/vatom.f \
	../m613/vlgndr.f \
	../m613/vrad.f \
	../m613/vylm.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m613

$(BINDIR)/m613: $(SRCS:.f=.o) $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m613

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m613.f
	cat $(SRCS) > ../source/m613.f

print: link
	lpr ../source/m613.f

clean:
	rm -f *.o core
