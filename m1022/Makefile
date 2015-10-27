#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1022
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
	../m1022/hfd2e.f \
	../m1022/m1022.f \
	../m1022/pm1022.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1022

$(BINDIR)/m1022: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1022

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1022.f
	cat $(SRCS) > ../source/m1022.f

print: link
	lpr ../source/m1022.f

clean:
	rm -f *.o core







