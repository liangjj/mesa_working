#   @(#)Makefile	1.2   11/28/95
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
TCGMSGLIB = /usr4/russo/src/pnl_stuff/tcgmsg/libtcgmsg.$(OPSYS).a
GET = sccs get
REL = 

SRCS = \
	../xm515/bakwrd.f \
	../xm515/bfgrd.f \
	../xm515/bigdij.f \
	../xm515/cfrt1.f \
	../xm515/cotes.f \
	../xm515/diis.f \
	../xm515/directj.f \
	../xm515/errmat.f \
	../xm515/estmq.f \
	../xm515/fmkclos.f \
	../xm515/fmkopen.f \
	../xm515/forwrd.f \
	../xm515/ftogrid.f \
	../xm515/ftogrid2.f \
	../xm515/ftoylm.f \
	../xm515/gdmat.f \
	../xm515/get1dm.f \
	../xm515/getgrid.f \
	../xm515/gridden.f \
	../xm515/gofish.f \
	../xm515/hunt.f \
	../xm515/jmatrix.f \
	../xm515/kmatrix.f \
	../xm515/kmtxguts.f \
	../xm515/kmopen.f \
	../xm515/ksenrg.f \
	../xm515/ksham.f \
	../xm515/laplace.f \
	../xm515/ldexp.f \
	../xm515/lebedev.f \
	../xm515/levshft.f \
	../xm515/m515.f \
	../xm515/mkatmg.f \
	../xm515/moden.f \
	../xm515/ncwts.f \
	../xm515/necote.f \
	../xm515/newton.f \
	../xm515/nrgwden.f \
	../xm515/nxtask.f \
	../xm515/pbasis.f \
	../xm515/plnkerr.f \
	../xm515/pm515.f \
	../xm515/poisson.f \
	../xm515/prefac.f \
	../xm515/prims.f \
	../xm515/prntld.f \
	../xm515/pseudo.f \
	../xm515/putj.f \
	../xm515/putvec.f \
	../xm515/qmat.f \
	../xm515/rsolver.f \
	../xm515/rhomax.f \
	../xm515/rules.f \
	../xm515/scttr.f \
	../xm515/setup.f \
	../xm515/sg1.f \
	../xm515/sizer.f \
	../xm515/sphere.f \
	../xm515/spline3.f \
	../xm515/splint.f \
	../xm515/sqzclo.f \
	../xm515/sqzopn.f \
	../xm515/symrsp.f \
	../xm515/tomo.f \
	../xm515/tr1dm.f \
	../xm515/tstylm.f \
	../xm515/v0.f \
	../xm515/v0atr.f \
	../xm515/vatom.f \
	../xm515/vylm.f \
	../xm515/vlgndr.f \
	../xm515/vrad.f \
	../xm515/vfmj.f \
	../xm515/vfmk.f \
	../xm515/vmkghi.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m515 $(BINDIR)/m515.slave

../xm515/getscm.small.o: ../xm515/getscm.small.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) ../xm515/getscm.small.$(OPSYS).f

../xm515/getscm.big.o: ../xm515/getscm.big.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) ../xm515/getscm.big.$(OPSYS).f

../xm515/getscm.tiny.o: ../xm515/getscm.tiny.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) ../xm515/getscm.tiny.$(OPSYS).f

$(BINDIR)/m515: $(SRCS:.f=.o) ../xm515/getscm.big.o $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../xm515/getscm.big.o $(MESALIB)  $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m515

$(BINDIR)/m515.slave: $(SRCS:.f=.o) ../xm515/getscm.small.o $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../xm515/getscm.small.o $(MESALIB)  $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m515.slave


$(BINDIR)/m515.slavelet: $(SRCS:.f=.o) ../xm515/getscm.tiny.o $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../xm515/getscm.tiny.o $(MESALIB)  $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m515.slavelet

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/xm515.f
	cat $(SRCS) > ../source/xm515.f

print: link
	lpr ../source/xm515.f

clean:
	rm -f *.o core


