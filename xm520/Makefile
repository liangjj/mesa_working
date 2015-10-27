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
	../m520/bakwrd.f \
	../m520/bfgrd.f \
	../m520/bigdij.f \
	../m520/cfrt1.f \
	../m520/cotes.f \
	../m520/diis.f \
	../m520/directj.f \
	../m520/errmat.f \
	../m520/estmq.f \
	../m520/fmkclos.f \
	../m520/fmkopen.f \
	../m520/forwrd.f \
	../m520/ftogrid.f \
	../m520/ftogrid2.f \
	../m520/ftoylm.f \
	../m520/getgrid.f \
	../m520/gdmat.f \
	../m520/get1dm.f \
	../m520/gridden.f \
	../m520/gofish.f \
	../m520/hunt.f \
	../m520/jmatrix.f \
	../m520/kmatrix.f \
	../m520/kmtxguts.f \
	../m520/kmopen.f \
	../m520/ksenrg.f \
	../m520/ksham.f \
	../m520/laplace.f \
	../m520/ldexp.f \
	../m520/levshft.f \
	../m520/m520.f \
	../m520/moden.f \
	../m520/ncwts.f \
	../m520/necote.f \
	../m520/newton.f \
	../m520/nrgwden.f \
	../m520/pm511.f \
	../m520/poisson.f \
	../m520/prefac.f \
	../m520/prims.f \
	../m520/prntld.f \
	../m520/pseudo.f \
	../m520/putj.f \
	../m520/putvec.f \
	../m520/qmat.f \
	../m520/rsolver.f \
	../m520/rules.f \
	../m520/scttr.f \
	../m520/setup.f \
	../m520/sizer.f \
	../m520/spline3.f \
	../m520/splint.f \
	../m520/sqzclo.f \
	../m520/sqzopn.f \
	../m520/symrsp.f \
	../m520/tomo.f \
	../m520/tr1dm.f \
	../m520/tstylm.f \
	../m520/v0.f \
	../m520/v0atr.f \
	../m520/vatom.f \
	../m520/vylm.f \
	../m520/vlgndr.f \
	../m520/vrad.f \
	../m520/vfmj.f \
	../m520/vmkghi.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m520

$(BINDIR)/m520: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB)  $(MDLIB) -o $(BINDIR)/m520

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m520.f
	cat $(SRCS) > ../source/m520.f

print: link
	lpr ../source/m520.f

clean:
	rm -f *.o core


