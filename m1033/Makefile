#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1033
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
	../m1033/cider.f \
	../m1033/m1033.f \
	../m1033/pm1033.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1033

$(BINDIR)/m1033: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1033

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1033.f
	cat $(SRCS) > ../source/m1033.f

print: link
	lpr ../source/m1033.f

clean:
	rm -f *.o core







