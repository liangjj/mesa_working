# @(#)Makefile	5.2 2/5/95
#
# Makefile for M721
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
	../m721/d1xc.f \
	../m721/dxcint.f \
	../m721/dxclos.f \
	../m721/dxopen.f \
	../m721/fmaccum.f \
	../m721/m721.f \
	../m721/pm721.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m721

$(BINDIR)/m721: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m721

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m721.f
	cat $(SRCS) > ../source/m721.f

print: link
	lpr ../source/m721.f

clean:
	rm -f *.o core


