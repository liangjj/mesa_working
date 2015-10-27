#   %W%   %G%
#
# Makefile for udipc
#
FC = g77
FFLAGS = -c -O -finit-local-zero -finit-local-zero

LD = g77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
BLASLIB= /System/Library/Frameworks/vecLib.framework/versions/A/vecLib
GET = sccs get
REL = 

SRCS = \
	../dipc/udipc.f \
	../dipc/dipint.f \
	../dipc/bforthog.f \
	../dipc/trnmat.f \
	../dipc/rdbinsqr.f \
	../dipc/rdbintri.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/xdipc

$(BINDIR)/xdipc: $(SRCS:.f=.o) $(MESALIB) $(BLASLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(BLASLIB) -o $(BINDIR)/xdipc

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/udipc.f
	cat $(SRCS) > ../source/udipc.f

print: link
	lpr ../source/udipc.f

clean:
	rm -f *.o core
