#    @(#)Makefile	5.1   11/6/94
#
#    Makefile for the symmetry library of MesaLib
#    
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR = $(UTILITY)
include $(TOP_DIR)/Makefile.inc
GET = sccs get
REL = 

SRCS = \
	atmset.f \
	bftoso.f \
	c2prim.f \
	ci.f \
	cirset.f \
	cnaxis.f \
	deornt.f \
	dmprep.f \
	equiv.f \
	findc2.f \
	findcn.f \
	findv.f \
	genrtr.f \
	group.f \
	grpsym.f \
	invert.f \
	move.f \
	octa.f \
	oe.f \
	oper.f \
	or3mom.f \
	oraxis.f \
	orc2v.f \
	orcn.f \
	ord2h.f \
	ordn.f \
	ordoc.f \
	orkey.f \
	ornax.f \
	orplan.f \
	orptst.f \
	oryz.f \
	ot2.f \
	permut.f \
	projct.f \
	ptgrp.f \
	putsym.f \
	qcentr.f \
	reflct.f \
	relate.f \
	rotate.f \
	rotx.f \
	roty.f \
	rotz.f \
	salc.f \
	secmom.f \
	sigmad.f \
	sigmah.f \
	sigmav.f \
	snaxis.f \
	sphere.f \
	sphset.f \
	sprint.f \
	symm.f \
	sympt.f \
	symset.f \
	tform.f \
	tmat.f \
	tolocl.f \
	tomast.f \
	triang.f \
	trmat.f \
	tstc3.f \
	tstc4.f \
	tstc5.f \
	vert.f 

OBJECTS = \
	$(SRCS:.f=.o) 

all: $(SYMLIB) $(OBJECTS)
$(SYMLIB): $(OBJECTS) 
	$(AR) $(ARFLAGS) $(SYMLIB) $(OBJECTS)
	$(RANLIB) $(SYMLIB)
	$(AR) $(ARFLAGS) $(UTILITY_LIB) $(OBJECTS)
	$(RANLIB) $(UTILITY_LIB)

sources: $(SRCS) 
$(SRCS):
	$(GET) $(RELEASE) $@

clean:
	rm -f *.o
	rm -f core
