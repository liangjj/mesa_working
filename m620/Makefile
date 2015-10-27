#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M620
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
	../m620/aintch.f \
	../m620/bond.f \
	../m620/center.f \
	../m620/clij.f \
	../m620/deriv.f \
	../m620/fdist.f \
	../m620/lonep.f \
	../m620/m620.f \
	../m620/norm.f \
	../m620/outrho.f \
	../m620/pdm88.f \
	../m620/pm620.f \
	../m620/sitout.f \
	../m620/syminv.f \
	../m620/symlin.f \
	../m620/vary.f \
	../m620/vcalc.f \
	../m620/vcalcr.f \
	../m620/vs.f \
	../m620/vxv.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m620

$(BINDIR)/m620: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m620

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m620.f
	cat $(SRCS) > ../source/m620.f

print: link
	lpr ../source/m620.f

clean:
	rm -f *.o core


