#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M732
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
	../m732/cztran.f \
	../m732/m732.f \
	../m732/pm732.f \
	../m732/trmain.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m732

$(BINDIR)/m732: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m732

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m732.f
	cat $(SRCS) > ../source/m732.f

print: link
	lpr ../source/m732.f

clean:
	rm -f *.o core


