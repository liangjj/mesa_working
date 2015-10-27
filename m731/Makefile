#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M731
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
	../m731/bend.f \
	../m731/fcorpr.f \
	../m731/formbg.f \
	../m731/fzprnt.f \
	../m731/m731.f \
	../m731/pm731.f \
	../m731/putf.f \
	../m731/rotf.f \
	../m731/rotff.f \
	../m731/rotfff.f \
	../m731/str.f \
	../m731/tors.f \
	../m731/tranf.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m731

$(BINDIR)/m731: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m731

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m731.f
	cat $(SRCS) > ../source/m731.f

print: link
	lpr ../source/m731.f

clean:
	rm -f *.o core


