#    @(#)Makefile	5.1   11/6/94
#
#    Makefile for the utility portion of MesaLib
#    
MESALIB = ../mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
GET = sccs get
REL = 

SRCS = \
	../util/basis.f \
	../util/chainx.f \
	../util/chkpnt.f \
	../util/corpr1.f \
	../util/dismat.f \
	../util/dollar.f \
	../util/drum.f \
	../util/enforc.f \
	../util/filchg.f \
	../util/fillel.f \
	../util/filmas.f \
	../util/fixdrt.f \
	../util/fmstrt.f \
	../util/frqprt.f \
	../util/funcnm.f \
	../util/getcf.f \
	../util/liopen.f \
	../util/masswt.f \
	../util/matout.f \
	../util/matprt.f \
	../util/norm1.f \
	../util/normal.f \
	../util/print.f \
	../util/project.f \
	../util/prsrte.f \
	../util/rcoord.f \
	../util/replic.f \
	../util/rtelbl.f \
	../util/rzmat.f \
	../util/rzsub.f \
	../util/subvar.f \
	../util/togenc.f \
	../util/vec.f \
	../util/vecout.f \
	../util/vprod.f \
	../util/wcoord.f \
	../util/wlmat.f \
	../util/wmat.f \
	../util/wvec.f \
	../util/wzmat.f \
	../util/wzsub.f \
	../util/zconv.f \
	../util/zprint.f \
	../util/ztoc.f

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
	rm -f ../../source/util.f
	cat $(SRCS) > ../../source/util.f

print: link
	lpr ../../source/util.f

clean:
	rm -f *.o
	rm -f core
