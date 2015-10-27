#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1951
#
FC = $FC)
FFLAGS = $(FFLAGS)

LD = $(FC)
LDFLAGS =

BINDIR = $(MESA_BIN)
MESALIB = $(MESA_LIB)/mesalib.a
MDLIB = $(MDLIB)
GET = sccs get
REL = 

SRCS = \
	m1951.f \
	mulpop.f \
	prpone.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1951

$(BINDIR)/m1951: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1951

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1951.f
	cat $(SRCS) > ../source/m1951.f

print: link
	lpr ../source/m1951.f

clean:
	rm -f *.o core


