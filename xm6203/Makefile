#  @(#)Makefile	1.2 10/27/94
#
# Makefile for M6203
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
	../m6203/banmat.f \
	../m6203/bkwrd.f \
	../m6203/bsplvb.f \
	../m6203/bsplpp.f \
	../m6203/cardin.f \
	../m6203/caspln.f \
	../m6203/convr.f \
	../m6203/copy.f \
	../m6203/drvint.f \
	../m6203/fact.f \
	../m6203/fbgint.f \
	../m6203/fbnint.f \
	../m6203/fndbrk.f \
	../m6203/forwrd.f \
	../m6203/gamfun.f \
	../m6203/interv.f \
	../m6203/legend.f \
	../m6203/legint.f \
	../m6203/lmdcmp.f \
	../m6203/m6203.f \
	../m6203/miscfn.f \
	../m6203/mkspln.f \
	../m6203/mkylm.f \
	../m6203/mvmul.f \
	../m6203/nrmylm.f \
	../m6203/plmtof.f \
	../m6203/pm6203.f \
	../m6203/prespl.f \
	../m6203/ppval.f \
	../m6203/prntfm.f \
	../m6203/rcbes.f \
	../m6203/runiqe.f \
	../m6203/scalyu.f \
	../m6203/scmphi.f \
	../m6203/splmat.f \
	../m6203/splcof.f \
	../m6203/sqrtwt.f \
	../m6203/srtleg.f \
	../m6203/tstzro.f \
	../m6203/vmmul.f \
	../m6203/ylmint.f \
	../m6203/yukawa.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6203

$(BINDIR)/m6203: $(SRCS:.f=.o) $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m6203

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m6203.f
	cat $(SRCS) > ../source/m6203.f

print: link
	lpr ../source/m6203.f

clean:
	rm -f *.o core
