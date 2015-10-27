#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1003
#
include $(MESA_HOME)/Makefile.inc

OTHERLIBS = \
	../library/m820lib.a \
	../library/m830lib.a \
	../library/m902lib.a \
	../library/m903lib.a
GET = sccs get
REL = 

SRCS = \
	../m1003/addlag.f\
	../m1003/ainex.f \
	../m1003/aint.f \
	../m1003/beapat.f\
	../m1003/blfold.f\
	../m1003/cint.f \
	../m1003/cpledr.f\
	../m1003/dcan.f\
	../m1003/derint.f\
	../m1003/derlag.f\
	../m1003/fixlag.f\
	../m1003/fixta.f\
	../m1003/force.f\
	../m1003/gethes.f\
	../m1003/getrab.f\
	../m1003/gttaao.f\
	../m1003/gttamo.f\
	../m1003/hesvec.f\
	../m1003/linear.f\
	../m1003/m1003.f\
	../m1003/maketa.f\
	../m1003/mccg1.f\
	../m1003/mccg1e.f\
	../m1003/mccga.f\
	../m1003/mccgb.f\
	../m1003/mccgc.f\
	../m1003/mccgpk.f\
	../m1003/mcfmtr.f\
	../m1003/mcg1.f\
	../m1003/mcg1e.f\
	../m1003/mcgcor.f\
	../m1003/mcgrd.f\
	../m1003/mcgupk.f\
	../m1003/mcgvca.f\
	../m1003/mcgvcb.f\
	../m1003/mcjcbi.f\
	../m1003/mcmxvc.f \
	../m1003/pm1003.f\
	../m1003/projct.f \
	../m1003/rdiagh.f\
	../m1003/scdm1.f\
	../m1003/setta.f\
	../m1003/tag1.f\
	../m1003/tag1e.f\
	../m1003/tagcor.f\
	../m1003/taint.f\
	../m1003/umat.f\
	../m1003/vsamul.f \
	../m1003/wtdiag.f \
	../m1003/xtalag.f\
	../m1003/yfold.f\
	../m1003/ytalag.f

all: $(MESA_BIN) $(SRCS:.f=.o)
	$(LOAD)/m1003

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1003.f
	cat $(SRCS) > ../source/m1003.f

print: link
	lpr ../source/m1003.f

clean:
	rm -f *.o core



