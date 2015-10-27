#   @(#)Makefile	4.1   7/7/93
#
# Makefile for M1990
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
	../m1990/dirac.f \
	../m1990/grid.f \
	../m1990/m1990.f \
	../m1990/pm1990.f \
	../m1990/xyzout.f 
 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1990
$(BINDIR)/m1990: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1990

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1990.f
	cat $(SRCS) > ../source/m1990.f

print: link
	lpr ../source/m1990.f

clean:
	rm -f *.o core



