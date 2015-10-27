#   %W%   %G%
#
# Makefile for M1001
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
OTHERLIBS = \
	../library/m820lib.a \
	../library/m830lib.a \
	../library/m902lib.a \
	../library/m903lib.a
GET = sccs get
REL = 

SRCS = \
	../m1001/addlag.f\
	../m1001/beapat.f\
	../m1001/blfold.f\
	../m1001/cpledr.f\
	../m1001/dcan.f\
	../m1001/derint.f\
	../m1001/derlag.f\
	../m1001/fixlag.f\
	../m1001/fixta.f\
	../m1001/force.f\
	../m1001/gethes.f\
	../m1001/getrab.f\
	../m1001/gttaao.f\
	../m1001/gttamo.f\
	../m1001/hesvec.f\
	../m1001/linear.f\
	../m1001/m1001.f\
	../m1001/maketa.f\
	../m1001/mccg1.f\
	../m1001/mccg1e.f\
	../m1001/mccga.f\
	../m1001/mccgb.f\
	../m1001/mccgc.f\
	../m1001/mccgpk.f\
	../m1001/mcfmtr.f\
	../m1001/mcg1.f\
	../m1001/mcg1e.f\
	../m1001/mcgcor.f\
	../m1001/mcgrd.f\
	../m1001/mcgupk.f\
	../m1001/mcgvca.f\
	../m1001/mcgvcb.f\
	../m1001/mcjcbi.f\
	../m1001/mcmxvc.f \
	../m1001/pm1001.f\
	../m1001/rdiagh.f\
	../m1001/scdm1.f\
	../m1001/setta.f\
	../m1001/tag1.f\
	../m1001/tag1e.f\
	../m1001/tagcor.f\
	../m1001/taint.f\
	../m1001/umat.f\
	../m1001/xtalag.f\
	../m1001/yfold.f\
	../m1001/ytalag.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1001

$(BINDIR)/m1001: $(SRCS:.f=.o) $(MESALIB) $(OTHERLIBS)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(OTHERLIBS) $(MESALIB)\
             -o $(BINDIR)/m1001

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1001.f
	cat $(SRCS) > ../source/m1001.f

print: link
	lpr ../source/m1001.f

clean:
	rm -f *.o core



