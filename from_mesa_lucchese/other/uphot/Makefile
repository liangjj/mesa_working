#   %W%   %G%
#
# Makefile for uphot
#
FC = g77
FFLAGS = -c -O -finit-local-zero -finit-local-zero

LD = g77
LDFLAGS =

BINDIR = ../bin
#MESALIB = ../library/mesalib.a
#BLASLIB= /System/Library/Frameworks/vecLib.framework/versions/A/vecLib
GET = sccs get
REL = 

SRCS = \
	../uphot/uphot.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/xphot

$(BINDIR)/xphot: $(SRCS:.f=.o) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  -o $(BINDIR)/xphot

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/uphot.f
	cat $(SRCS) > ../source/uphot.f

print: link
	lpr ../source/uphot.f

clean:
	rm -f *.o core
