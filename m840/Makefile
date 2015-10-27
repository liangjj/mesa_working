#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M840
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
	../m840/m840.f \
	../m840/mn840.f \
	../m840/pm840.f \
	../m840/toder.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m840

$(BINDIR)/m840: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m840

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m840.f
	cat $(SRCS) > ../source/m840.f

print: link
	lpr ../source/m840.f

clean:
	rm -f *.o core




