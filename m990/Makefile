#   @(#)Makefile	1.4   8/1/91
#
# Makefile for M990
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../m990/kzero.f \
	../m990/m990.f \
	../m990/pm990.f \
	../m990/prwalk.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m990

$(BINDIR)/m990: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/m990

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m990.f
	cat $(SRCS) > ../source/m990.f

print: link
	lpr ../source/m990.f

clean:
	rm -f *.o core


