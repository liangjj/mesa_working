#   %W%   %G%
#
# Makefile for pffbfg
#
FC = g77
FFLAGS = -c -O -finit-local-zero -finit-local-zero

LD = g77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
BLASLIB = /usr/lib/libblas.a
GET = sccs get
REL = 

SRCS = \
	../pffbfg/bfints.f \
	../pffbfg/pffbfg.f \
	../pffbfg/ffints.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/xpffbfg

$(BINDIR)/xpffbfg: $(SRCS:.f=.o) $(MESALIB) $(BLASLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(BLASLIB) -o $(BINDIR)/xpffbfg

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/cffbfg.f
	cat $(SRCS) > ../source/cffbfg.f

print: link
	lpr ../source/cffbfg.f

clean:
	rm -f *.o core
