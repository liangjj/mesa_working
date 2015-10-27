#   @(#)Makefile	1.1   11/30/90
#
# Makefile for M950
#
FC = f77
FFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../m950/basout.f \
	../m950/m950.f \
	../m950/pm950.f \
	../m950/wrbinsqr.f \
	../m950/wrbintri.f

.f.o:
	$(FC) $(FFLAGS) -c $<

all: $(BINDIR)/m950

$(BINDIR)/m950: $(SRCS:.f=.o) $(MESALIB)
	$(FC) $(FFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/m950

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m950.f
	cat $(SRCS) > ../source/m950.f

print: link
	lpr ../source/m950.f

clean:
	rm -f *.o core



