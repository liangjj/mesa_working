#   @(#)Makefile	5.2   11/28/95
#
# Makefile for M515
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB =
#  temporary fix for the sgi
#TCGMSGLIB = $(MESA_LIBS)/tcgmsg/libtcgmsg.$(OPSYS).a
TCGMSGLIB = ../library/tcgmsg/libtcgmsg.$(OPSYS).a
GET = sccs get
REL = 

SRCS = \
	../m515/bakwrd.f \
	../m515/bfgrd.f \
	../m515/bigdij.f \
	../m515/cfrt1.f \
	../m515/cotes.f \
	../m515/diis.f \
	../m515/directj.f \
	../m515/errmat.f \
	../m515/estmq.f \
	../m515/fmkclos.f \
	../m515/fmkopen.f \
	../m515/forwrd.f \
	../m515/ftogrid.f \
	../m515/ftogrid2.f \
	../m515/ftoylm.f \
	../m515/gdmat.f \
	../m515/get1dm.f \
	../m515/getgrid.f \
	../m515/gridden.f \
	../m515/gofish.f \
	../m515/hunt.f \
	../m515/jmatrix.f \
	../m515/kmatrix.f \
	../m515/kmtxguts.f \
	../m515/kmopen.f \
	../m515/ksenrg.f \
	../m515/ksham.f \
	../m515/laplace.f \
	../m515/ldexp.f \
	../m515/lebedev.f \
	../m515/levshft.f \
	../m515/m515.f \
	../m515/mkatmg.f \
	../m515/moden.f \
	../m515/ncwts.f \
	../m515/necote.f \
	../m515/newton.f \
	../m515/nrgwden.f \
	../m515/nxtask.f \
	../m515/pbasis.f \
	../m515/plnkerr.f \
	../m515/pm515.f \
	../m515/poisson.f \
	../m515/prefac.f \
	../m515/prims.f \
	../m515/prntld.f \
	../m515/pseudo.f \
	../m515/putj.f \
	../m515/putvec.f \
	../m515/qmat.f \
	../m515/rsolver.f \
	../m515/rhomax.f \
	../m515/rules.f \
	../m515/scttr.f \
	../m515/setup.f \
	../m515/sg1.f \
	../m515/sizer.f \
	../m515/sphere.f \
	../m515/spline3.f \
	../m515/splint.f \
	../m515/sqzclo.f \
	../m515/sqzopn.f \
	../m515/symrsp.f \
	../m515/tomo.f \
	../m515/tr1dm.f \
	../m515/tstylm.f \
	../m515/v0.f \
	../m515/v0atr.f \
	../m515/vatom.f \
	../m515/vylm.f \
	../m515/vlgndr.f \
	../m515/vrad.f \
	../m515/vfmj.f \
	../m515/vfmk.f \
	../m515/vmkghi.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m515 $(BINDIR)/m515.slave

../m515/getscm.small.o: ../m515/getscm.small.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) ../m515/getscm.small.$(OPSYS).f

../m515/getscm.big.o: ../m515/getscm.big.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) ../m515/getscm.big.$(OPSYS).f

../m515/getscm.tiny.o: ../m515/getscm.tiny.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) ../m515/getscm.tiny.$(OPSYS).f

$(BINDIR)/m515: $(SRCS:.f=.o) ../m515/getscm.big.o $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../m515/getscm.big.o $(MESALIB)  $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m515

$(BINDIR)/m515.slave: $(SRCS:.f=.o) ../m515/getscm.small.o $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../m515/getscm.small.o $(MESALIB)  $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m515.slave


$(BINDIR)/m515.slavelet: $(SRCS:.f=.o) ../m515/getscm.tiny.o $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../m515/getscm.tiny.o $(MESALIB)  $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m515.slavelet

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m515.f
	cat $(SRCS) > ../source/m515.f

print: link
	lpr ../source/m515.f

clean:
	rm -f *.o core


