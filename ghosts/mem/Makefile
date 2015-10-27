#   %W% %G%
#
#    Makefile for the memory manager library of MesaLib
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
	../memory/getpt.f \
	../memory/inimem.f \
	../memory/memdat.f \
	../memory/relpt.f
 
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
	rm -f ../../source/memory.f
	cat $(SRCS) > ../../source/memory.f

print: link
	lpr ../../source/memory.f

clean:
	rm -f *.o
	rm -f core
