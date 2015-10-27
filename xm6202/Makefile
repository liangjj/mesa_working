#  %W% %G%
#
# Makefile for M6202
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
	../m6202/bkwrd.f \
	../m6202/caspln.f \
	../m6202/convr.f \
	../m6202/copy.f \
	../m6202/drvint.f \
	../m6202/fact.f \
	../m6202/fbgint.f \
	../m6202/fbnint.f \
	../m6202/fndbrk.f \
	../m6202/forwrd.f \
	../m6202/homog.f \
	../m6202/legend.f \
	../m6202/legint.f \
	../m6202/lmdcmp.f \
	../m6202/miscfn.f \
	../m6202/mkspln.f \
	../m6202/mkylm.f \
	../m6202/m6202.f \
	../m6202/nrmylm.f \
	../m6202/plmtof.f \
	../m6202/pm6202.f \
	../m6202/prntfm.f \
	../m6202/runiqe.f \
	../m6202/scalyu.f \
	../m6202/scmphi.f \
	../m6202/setzro.f \
	../m6202/spline.f \
	../m6202/sqrtwt.f \
	../m6202/srtleg.f \
	../m6202/vmmul.f \
	../m6202/ylmint.f \
	../m6202/yukawa.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6202

$(BINDIR)/m6202: $(SRCS:.f=.o) $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m6202

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m6202.f
	cat $(SRCS) > ../source/m6202.f

print: link
	lpr ../source/m6202.f

clean:
	rm -f *.o core
