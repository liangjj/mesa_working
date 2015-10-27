#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M618
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
	../m618/getrad.f \
	../m618/m618.f \
	../m618/pm618.f \
	../m618/sobseq1.f \
	../m618/surface.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m618

$(BINDIR)/m618: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m618

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m618.f
	cat $(SRCS) > ../source/m618.f

print: link
	lpr ../source/m618.f

clean:
	rm -f *.o core


