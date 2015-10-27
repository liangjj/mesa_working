#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M701
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
	../m701/m701.f \
	../m701/pm701.f \
	../m701/tr1dm.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m701

$(BINDIR)/m701: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m701

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m701.f
	cat $(SRCS) > ../source/m701.f

print: link
	lpr ../source/m701.f

clean:
	rm -f *.o core



