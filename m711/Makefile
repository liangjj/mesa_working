#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M711
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
	../m711/d2expn.f \
	../m711/derexp.f \
	../m711/dprim.f \
	../m711/driver.f \
	../m711/fd2ant.f \
	../m711/fd2bnt.f \
	../m711/fmdint.f \
	../m711/get1dm.f \
	../m711/getalp.f \
	../m711/gj2pdm.f \
	../m711/ldexp.f \
	../m711/m711.f \
	../m711/makd2g.f \
	../m711/makd2h.f \
	../m711/makd2i.f \
	../m711/makdg.f \
	../m711/makdh.f \
	../m711/makdi.f \
	../m711/md2ghi.f \
	../m711/mkdghi.f \
	../m711/movd2e.f \
	../m711/movder.f \
	../m711/pm711.f \
	../m711/prefac.f \
	../m711/prims.f \
	../m711/prntld.f \
	../m711/redund.f \
	../m711/sizer.f \
	../m711/vfmint.f \
	../m711/vmkghi.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m711

$(BINDIR)/m711: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m711

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m711.f
	cat $(SRCS) > ../source/m711.f

print: link
	lpr ../source/m711.f

clean:
	rm -f *.o core


