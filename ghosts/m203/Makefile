#   %W%   %G%
#
# Makefile for M203
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
	../m203/cztran.f \
	../m203/m203.f \
	../m203/pm203.f \
	../m203/trmain.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m203

$(BINDIR)/m203: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m203

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m203.f
	cat $(SRCS) > ../source/m203.f

print: link
	lpr ../source/m203.f

clean:
	rm -f *.o core


