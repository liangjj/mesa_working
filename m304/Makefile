#      @(#)Makefile	5.1   11/6/94
#
#      Makefile for M304
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
	../m304/m304.f \
	../m304/soint.f \
	../m304/pm304.f \
	../m304/put1so.f \
	../m304/sox.f \
	../m304/soy.f \
	../m304/soz.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m304

$(BINDIR)/m304: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m304

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m304.f
	cat $(SRCS) > ../source/m304.f

print: link
	lpr ../source/m304.f

clean:
	rm -f *.o core


