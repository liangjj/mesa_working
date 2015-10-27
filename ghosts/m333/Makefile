#   @(#)Makefile	1.1   11/30/90
#
# Makefile for M333
#
FC = f77
FFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../m333/m333.f \
	../m333/mn333.f \
	../m333/pm333.f \
	../m333/sort32.f \
	../m333/sort64.f 

.f.o:
	$(FC) $(FFLAGS) -c $<

all: $(BINDIR)/m333

$(BINDIR)/m333: $(SRCS:.f=.o) $(MESALIB)
	$(FC) $(FFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/m333

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m333.f
	cat $(SRCS) > ../source/m333.f

print: link
	lpr ../source/m333.f

clean:
	rm -f *.o core




