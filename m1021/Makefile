#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1021
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
	../m1021/dumbu.f \
	../m1021/fillu.f \
	../m1021/m1021.f \
	../m1021/mcscfu.f \
	../m1021/pm1021.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1021

$(BINDIR)/m1021: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1021

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1021.f
	cat $(SRCS) > ../source/m1021.f

print: link
	lpr ../source/m1021.f

clean:
	rm -f *.o core







