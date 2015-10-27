#    @(#)Makefile	5.2   4/17/95
#
#    Makefile for the mathematics library of MesaLib
#    
MESALIB = ../mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

GET = sccs get
REL = 

SRCS = \
	../math/ambc.f \
	../math/ambct.f \
	../math/ambtc.f \
	../math/ambtct.f \
	../math/apbc.f \
	../math/apbct.f \
	../math/apbctt.f \
	../math/apbtc.f \
	../math/apbtct.f \
	../math/arrmax.f \
	../math/atebc.f \
	../math/atembc.f \
	../math/atmbc.f \
	../math/atpbc.f \
	../math/atpbct.f \
	../math/bico.f \
	../math/bliu.f \
	../math/degrsp.f \
	../math/ebc.f \
	../math/ebct.f \
	../math/ebtc.f \
	../math/ebtct.f \
	../math/embc.f \
	../math/embct.f \
	../math/embtc.f \
	../math/embtct.f \
	../math/factln.f \
	../math/factrl.f \
	../math/flin.f \
	../math/fold.f \
	../math/gammln.f \
	../math/givens.f \
	../math/gvadd.f \
	../math/gvdiv.f \
	../math/gvmul.f \
	../math/gvsub.f \
	../math/iveceq.f \
	../math/izero.f \
	../math/lowdin.f \
	../math/minvrt.f \
	../math/runit.f \
	../math/rveceq.f \
	../math/rzero.f \
	../math/sadd.f \
	../math/schmdt.f \
	../math/schmidt.f \
	../math/sgmm.f \
	../math/sinv.f \
	../math/smul.f \
	../math/solver.f \
	../math/sorter.f \
	../math/sorter1.f \
	../math/sorter2.f \
	../math/sorter3.f \
	../math/sorter4.f \
	../math/sqtotr.f \
	../math/ssub.f \
	../math/trace.f \
	../math/trtosq.f \
	../math/vadd.f \
	../math/vclean.f \
	../math/vdiv.f \
	../math/vexp.f \
	../math/vfill.f \
	../math/vinv.f \
	../math/vmove.f \
	../math/vmul.f \
	../math/vadmul.f \
	../math/vneg.f \
	../math/vpoly.f \
	../math/vpolyd.f \
	../math/vpower.f \
	../math/vsqrt.f \
	../math/vsub.f \
	../math/vwx.f \
	../math/vwxs.f \
	../math/vwxy.f \
	../math/vwxyi.f 

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
	rm -f ../../source/math.f
	cat $(SRCS) > ../../source/math.f

print: link
	lpr ../../source/math.f

clean:
	rm -f *.o
	rm -f core
