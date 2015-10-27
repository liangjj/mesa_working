#    @(#)Makefile	5.1   11/6/94
#
#    Makefile for the integral library of MesaLib
#    
#    
MESALIB = ../mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

GET = sccs get
REL = 

SRCS = \
	../intgrl/ang1.f \
	../intgrl/ang2.f \
	../intgrl/bess.f \
	../intgrl/bessv.f \
	../intgrl/cfrt1.f \
	../intgrl/cfrt2.f \
	../intgrl/derexp.f \
	../intgrl/facab.f \
	../intgrl/fmonel.f \
	../intgrl/fmt.f \
	../intgrl/kmprsz.f \
	../intgrl/ldata.f \
	../intgrl/ldexp.f \
	../intgrl/ltab.f \
	../intgrl/out32.f \
	../intgrl/out64.f \
	../intgrl/pseud1.f \
	../intgrl/pseud2.f \
	../intgrl/ptprep.f \
	../intgrl/put1el.f \
	../intgrl/qalt.f \
	../intgrl/qasy.f \
	../intgrl/qcomp.f \
	../intgrl/qpasy.f \
	../intgrl/qpow.f \
	../intgrl/qqrad.f \
	../intgrl/quadr.f \
	../intgrl/recur1.f \
	../intgrl/recur2.f \
	../intgrl/recurf.f \
	../intgrl/root4.f \
	../intgrl/root5.f \
	../intgrl/roots.f \
	../intgrl/rt123.f \
	../intgrl/sints.f \
	../intgrl/sizer1.f \
	../intgrl/sizer2.f \
	../intgrl/srchrv.f \
	../intgrl/stwod.f \
	../intgrl/trans.f \
	../intgrl/trans1.f \
	../intgrl/trans2.f \
	../intgrl/trssss.f \
	../intgrl/vfmint.f \
	../intgrl/vmakg.f \
	../intgrl/vmakh.f \
	../intgrl/vmaki2.f \
	../intgrl/vmkghi.f \
	../intgrl/vrt1.f \
	../intgrl/vrt2.f \
	../intgrl/xpandz.f \
	../intgrl/ztab.f
 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(MESALIB)
$(MESALIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(MESALIB) $(?:.f=.o)
	$(RANLIB) $(MESALIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@

link: sources
	rm -f ../../source/intgrl.f
	cat $(SRCS) > ../../source/intgrl.f

print: link
	lpr ../../source/intgrl.f

clean:
	rm -f *.o
	rm -f core
