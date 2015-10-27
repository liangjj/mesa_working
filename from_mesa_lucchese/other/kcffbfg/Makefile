#   %W%   %G%
#
# Makefile for cffbfg
#
# FC = g77
FFLAGS = $(FFLAGSOZ)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
BLASLIB= $(BLASUSE)
GET = sccs get
REL = 

SRCS = \
	../kcffbfg/bfints.f \
	../kcffbfg/cffbfg.f \
	../kcffbfg/ffints.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/xcffbfg

$(BINDIR)/xcffbfg: $(SRCS:.f=.o)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(BLASLIB) -o $(BINDIR)/xcffbfg

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
