#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1991
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
	../m1991/m1991.f \
	../m1991/mulpop.f \
	../m1991/pm1991.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1991

$(BINDIR)/m1991: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1991

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1991.f
	cat $(SRCS) > ../source/m1991.f

print: link
	lpr ../source/m1991.f

clean:
	rm -f *.o core


