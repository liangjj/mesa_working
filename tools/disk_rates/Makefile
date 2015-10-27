#   %W%   %G%
#
# Makefile for M1992
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
	m1992.f \
	pm1992.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1992

$(BINDIR)/m1992: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1992

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@
link: sources	
	rm -f ../source/m1992.f
	cat $(SRCS) > ../source/m1992.f
print: link
	lpr ../source/m1992.f
clean:
	rm -f *.o core


