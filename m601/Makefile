#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M601
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
	../m601/symprj.f \
	../m601/m601.f \
	../m601/pm601.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m601

$(BINDIR)/m601: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m601

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m601.f
	cat $(SRCS) > ../source/m601.f

print: link
	lpr ../source/m601.f

clean:
	rm -f *.o core


