#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M621
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
	../m621/bsj.f \
	../m621/input.f \
	../m621/lubksb.f \
	../m621/ludcmp.f \
	../m621/m621.f \
	../m621/pm621.f \
	../m621/setup.f \
	../m621/sobseq1.f \
	../m621/surface.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m621

$(BINDIR)/m621: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m621

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m621.f
	cat $(SRCS) > ../source/m621.f

print: link
	lpr ../source/m621.f

clean:
	rm -f *.o core


