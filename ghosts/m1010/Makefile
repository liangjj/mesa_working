#   @(#)Makefile	1.1   11/30/90
#
# Makefile for M1010
#
FC = f77
FFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../m1010/dfock.f \
	../m1010/lagrng.f \
	../m1010/m1010.f \
	../m1010/pm1010.f \
	../m1010/nucrep.f \
	../m1010/sder.f 

.f.o:
	$(FC) $(FFLAGS) -c $<

all: $(BINDIR)/m1010

$(BINDIR)/m1010: $(SRCS:.f=.o) $(MESALIB)
	$(FC) $(FFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/m1010

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1010.f
	cat $(SRCS) > ../source/m1010.f

print: link
	lpr ../source/m1010.f

clean:
	rm -f *.o core



