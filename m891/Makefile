#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M891
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
	../m891/pm891.f \
	../m891/m891.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m891

$(BINDIR)/m891: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m891

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m891.f
	cat $(SRCS) > ../source/m891.f

print: link
	lpr ../source/m891.f

clean:
	rm -f *.o core


