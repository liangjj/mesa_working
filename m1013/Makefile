#
#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1013
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
	../m1013/m1013.f \
	../m1013/pm1013.f \
	../m1013/makeb0.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1013

$(BINDIR)/m1013: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1013

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1013.f
	cat $(SRCS) > ../source/m1013.f

print: link
	lpr ../source/m1013.f

clean:
	rm -f *.o core






