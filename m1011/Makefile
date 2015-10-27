#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1011
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
	../m1011/jandks.f \
	../m1011/gdmat.f \
	../m1011/lagrng.f \
	../m1011/m1011.f \
	../m1011/pm1011.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1011

$(BINDIR)/m1011: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1011

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1011.f
	cat $(SRCS) > ../source/m1011.f

print: link
	lpr ../source/m1011.f

clean:
	rm -f *.o core




