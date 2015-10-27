#   @(#)Makefile	5.3 4/18/95
#
# Makefile for M511
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
	../m511/bakwrd.f \
	../m511/bigdij.f \
	../m511/cfrt1.f \
	../m511/cotes.f \
	../m511/diis.f \
	../m511/directj.f \
	../m511/errmat.f \
	../m511/estmq.f \
	../m511/fmkclos.f \
	../m511/fmkopen.f \
	../m511/forwrd.f \
	../m511/ftogrid.f \
	../m511/ftoylm.f \
	../m511/gdmat.f \
	../m511/get1dm.f \
	../m511/gofish.f \
	../m511/gofish2.f \
	../m511/interv.f \
	../m511/jmatrix.f \
	../m511/kmatrix.f \
	../m511/kmtxguts.f \
	../m511/kmopen.f \
	../m511/ksenrg.f \
	../m511/ksham.f \
	../m511/laplace.f \
	../m511/ldexp.f \
	../m511/levshft.f \
	../m511/m511.f \
	../m511/ncwts.f \
	../m511/necote.f \
	../m511/newton.f \
	../m511/nrgwden.f \
	../m511/pm511.f \
	../m511/poisson.f \
	../m511/prefac.f \
	../m511/prims.f \
	../m511/prntld.f \
	../m511/pseudo.f \
	../m511/putj.f \
	../m511/putvec.f \
	../m511/qmat.f \
	../m511/rsolver.f \
	../m511/rules.f \
	../m511/scttr.f \
	../m511/setup.f \
	../m511/sizer.f \
	../m511/spline3.f \
	../m511/splint.f \
	../m511/sqzclo.f \
	../m511/sqzopn.f \
	../m511/symrsp.f \
	../m511/tomo.f \
	../m511/tr1dm.f \
	../m511/tstylm.f \
	../m511/v0.f \
	../m511/v0atr.f \
	../m511/vatom.f \
	../m511/vylm.f \
	../m511/vlgndr.f \
	../m511/vrad.f \
	../m511/vfmj.f \
	../m511/vmkghi.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m511

$(BINDIR)/m511: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB)  $(MDLIB) -o $(BINDIR)/m511

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m511.f
	cat $(SRCS) > ../source/m511.f

print: link
	lpr ../source/m511.f

clean:
	rm -f *.o core


