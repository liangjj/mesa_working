#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M4001
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
	../m4001/cogen.f \
	../m4001/derfc.f \
	../m4001/m4001.f \
	../m4001/symm.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m4001

$(BINDIR)/m4001: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m4001

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m4001.f
	cat $(SRCS) > ../source/m4001.f

print: link
	lpr ../source/m4001.f

clean:
	rm -f *.o core
