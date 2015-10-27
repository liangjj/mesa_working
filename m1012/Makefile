#   @(#)Makefile	5.1   11/6/94
# Makefile for M1012
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
	../m1012/hessn.f \
	../m1012/m1012.f \
	../m1012/pm1012.f \
	../m1012/ptrs.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1012

$(BINDIR)/m1012: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1012

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1012.f
	cat $(SRCS) > ../source/m1012.f

print: link
	lpr ../source/m1012.f

clean:
	rm -f *.o core






