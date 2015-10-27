#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M619
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
	../m619/gcoord.f \
	../m619/m619.f \
	../m619/pm619.f \
	../m619/v0.f \
	../m619/v0atr.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m619

$(BINDIR)/m619: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m619

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m619.f
	cat $(SRCS) > ../source/m619.f

print: link
	lpr ../source/m619.f

clean:
	rm -f *.o core


