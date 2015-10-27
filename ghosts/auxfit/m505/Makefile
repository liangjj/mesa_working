#   %W%   %G%
#
# Makefile for M505
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../xm505/diis.f \
	../xm505/fixlag.f \
	../xm505/formjk.f \
	../xm505/gdmat.f \
	../xm505/gvbcof.f \
	../xm505/hfenrg.f \
	../xm505/khnvec.f \
	../xm505/lagrng.f \
	../xm505/m505.f \
	../xm505/pm505.f \
	../xm505/pseud.f \
	../xm505/putvec.f \
	../xm505/qmat.f \
	../xm505/scf.f \
	../xm505/setup.f \
	../xm505/sortbl.f \
	../xm505/srtv.f \
	../xm505/symrsp.f \
	../xm505/unique.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m505

$(BINDIR)/m505: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/m505

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m505.f
	cat $(SRCS) > ../source/m505.f

print: link
	lpr ../source/m505.f

clean:
	rm -f *.o core


