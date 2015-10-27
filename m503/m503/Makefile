#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M503
#
FC = fort
FFLAGS = -c

LD = fort
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB = /usr/lib/libcxml.a
GET = sccs get
REL = 

SRCS = \
	../m503/diis.f \
	../m503/fixlag.f \
	../m503/gdmat.f \
	../m503/gvbcof.f \
	../m503/hfenrg.f \
	../m503/jandks.f \
	../m503/khnvec.f \
	../m503/lagrng.f \
	../m503/m503.f \
	../m503/pm503.f \
	../m503/pseud.f \
	../m503/putvec.f \
	../m503/qmat.f \
	../m503/scf.f \
	../m503/setup.f \
	../m503/symrsp.f \
	../m503/tomo.f \
	../m503/tmptou.f \
	../m503/unique.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m503

$(BINDIR)/m503: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m503

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m503.f
	cat $(SRCS) > ../source/m503.f

print: link
	lpr ../source/m503.f

clean:
	rm -f *.o core


