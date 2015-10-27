#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1020
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
	../m1020/direct.f \
	../m1020/m1020.f \
	../m1020/pm1020.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1020

$(BINDIR)/m1020: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1020

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1020.f
	cat $(SRCS) > ../source/m1020.f

print: link
	lpr ../source/m1020.f

clean:
	rm -f *.o core







