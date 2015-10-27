#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M712
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB =
GET = sccs get
REL = 

SRCS = \
	../m712/ci2dm.f \
	../m712/d2expn.f \
	../m712/derexp.f \
	../m712/dprim.f \
	../m712/driver.f \
	../m712/dsort.f \
	../m712/fd2ant.f \
	../m712/fd2bnt.f \
	../m712/fmdint.f \
	../m712/get1dm.f \
	../m712/getalp.f \
	../m712/ghf2dm.f \
	../m712/hf2pdm.f \
	../m712/ldexp.f \
	../m712/m712.f \
	../m712/makd2g.f \
	../m712/makd2h.f \
	../m712/makd2i.f \
	../m712/makdg.f \
	../m712/makdh.f \
	../m712/makdi.f \
	../m712/md2ghi.f \
	../m712/mkdghi.f \
	../m712/movd2e.f \
	../m712/movder.f \
	../m712/nomgrp.f \
	../m712/pm712.f \
	../m712/prefac.f \
	../m712/prims.f \
	../m712/prntld.f \
	../m712/redund.f \
	../m712/sizer.f \
	../m712/thf2dm.f \
	../m712/tmc2dm.f \
	../m712/tr1dm.f \
	../m712/trtpdm.f \
	../m712/vfmint.f \
	../m712/vmkghi.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m712

$(BINDIR)/m712: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m712

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m712.f
	cat $(SRCS) > ../source/m712.f

print: link
	lpr ../source/m712.f

clean:
	rm -f *.o core


