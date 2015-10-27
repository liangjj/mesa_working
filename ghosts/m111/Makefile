#   %W%   %G%
#
# Makefile for M111
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
	../m111/cztran.f \
	../m111/m111.f \
	../m111/pm111.f \
	../m111/trmain.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m111

$(BINDIR)/m111: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m111

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m111.f
	cat $(SRCS) > ../source/m111.f

print: link
	lpr ../source/m111.f

clean:
	rm -f *.o core


