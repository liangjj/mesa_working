#    @(#)Makefile	5.1   11/6/94
#
#    Makefile for the symmetry library of MesaLib
#    
MESALIB = ../mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

GET = 
REL = 

SRCS = \
	../sym/atmset.f \
	../sym/bftoso.f \
	../sym/c2prim.f \
	../sym/ci.f \
	../sym/cirset.f \
	../sym/cnaxis.f \
	../sym/deornt.f \
	../sym/dmprep.f \
	../sym/equiv.f \
	../sym/findc2.f \
	../sym/findcn.f \
	../sym/findv.f \
	../sym/genrtr.f \
	../sym/group.f \
	../sym/grpsym.f \
	../sym/invert.f \
	../sym/move.f \
	../sym/octa.f \
	../sym/oe.f \
	../sym/oper.f \
	../sym/or3mom.f \
	../sym/oraxis.f \
	../sym/orc2v.f \
	../sym/orcn.f \
	../sym/ord2h.f \
	../sym/ordn.f \
	../sym/ordoc.f \
	../sym/orkey.f \
	../sym/ornax.f \
	../sym/orplan.f \
	../sym/orptst.f \
	../sym/oryz.f \
	../sym/ot2.f \
	../sym/permut.f \
	../sym/projct.f \
	../sym/ptgrp.f \
	../sym/putsym.f \
	../sym/qcentr.f \
	../sym/reflct.f \
	../sym/relate.f \
	../sym/rotate.f \
	../sym/rotx.f \
	../sym/roty.f \
	../sym/rotz.f \
	../sym/salc.f \
	../sym/secmom.f \
	../sym/sigmad.f \
	../sym/sigmah.f \
	../sym/sigmav.f \
	../sym/snaxis.f \
	../sym/sphere.f \
	../sym/sphset.f \
	../sym/sprint.f \
	../sym/symm.f \
	../sym/sympt.f \
	../sym/symset.f \
	../sym/tform.f \
	../sym/tmat.f \
	../sym/tolocl.f \
	../sym/tomast.f \
	../sym/triang.f \
	../sym/trmat.f \
	../sym/tstc3.f \
	../sym/tstc4.f \
	../sym/tstc5.f \
	../sym/vert.f 

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
	rm -f ../../source/sym.f
	cat $(SRCS) > ../../source/sym.f

print: link
	lpr ../../source/sym.f

clean:
	rm -f *.o
	rm -f core
