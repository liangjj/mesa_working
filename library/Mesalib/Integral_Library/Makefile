#    @(#)Makefile	5.1   11/6/94
#
#    Makefile for the integral library of MesaLib
#    
#    
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR = $(MESALIB)
include $(TOP_DIR)/Makefile.inc
GET = sccs get
REL = 
SRCS = \
	ang1.f \
	ang2.f \
	bess.f \
	bessv.f \
	cfrt1.f \
	cfrt2.f \
	derexp.f \
	facab.f \
	fmonel.f \
	fmt.f \
	kmprsz.f \
	ldata.f \
	ldexp.f \
	ltab.f \
	out32.f \
	out64.f \
	pseud1.f \
	pseud2.f \
	ptprep.f \
	put1el.f \
	qalt.f \
	qasy.f \
	qcomp.f \
	qpasy.f \
	qpow.f \
	qqrad.f \
	quadr.f \
	recur1.f \
	recur2.f \
	recurf.f \
	root4.f \
	root5.f \
	roots.f \
	rt123.f \
	sints.f \
	sizer1.f \
	sizer2.f \
	srchrv.f \
	stwod.f \
	trans.f \
	trans1.f \
	trans2.f \
	trssss.f \
	vfmint.f \
	vmakg.f \
	vmakh.f \
	vmaki2.f \
	vmkghi.f \
	vrt1.f \
	vrt2.f \
	xpandz.f \
	ztab.f
OBJECTS = \
	$(SRCS:.f=.o) 

all: $(INTLIB)
$(INTLIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $(INTLIB) $(OBJECTS)
	$(RANLIB) $(INTLIB)
	$(AR) $(ARFLAGS) $(MESA_LIB) $(OBJECTS)
	$(RANLIB) $(MESA_LIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@
clean:
	rm -f *.o
	rm -f core
