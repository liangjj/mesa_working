#   @(#)Makefile	5.2 2/5/95
#
#    Makefile for the density functional thoery library of MesaLib
#    
MESALIB = ../mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = f77
FFLAGS = -c
GET = sccs get
REL = 

SRCS = \
	../dft/alrixq.f \
	../dft/angsiz.f \
	../dft/atgrd.f \
	../dft/beckef.f \
	../dft/bfgrd.f \
	../dft/catenat.f \
	../dft/chebyq.f \
	../dft/eulerq.f \
	../dft/gthr.f \
	../dft/gauleg.f \
	../dft/gblk.f \
	../dft/gridden.f \
	../dft/lebedev.f \
	../dft/lgndrx.f \
	../dft/lypf.f \
	../dft/mkatmg.f \
	../dft/mkgrid.f \
	../dft/nzident.f \
	../dft/pakgrd.f \
	../dft/radial.f \
	../dft/rhomax.f \
	../dft/sg1.f \
	../dft/slaterf.f \
	../dft/sphere.f \
	../dft/tagrid.f \
        ../dft/vorgrad.f \
	../dft/voronoi.f \
	../dft/vwnf.f 
 
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
	rm -f ../../source/dft.f
	cat $(SRCS) > ../../source/dft.f

print: link
	lpr ../../source/dft.f

clean:
	rm -f *.o
	rm -f core
