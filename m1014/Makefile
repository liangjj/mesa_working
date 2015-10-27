#
#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1014
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
	../m1014/dumbwa.f \
	../m1014/m1014.f \
	../m1014/pm1014.f \
	../m1014/makewa.f 


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1014

$(BINDIR)/m1014: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1014

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1014.f
	cat $(SRCS) > ../source/m1014.f

print: link
	lpr ../source/m1014.f

clean:
	rm -f *.o core






