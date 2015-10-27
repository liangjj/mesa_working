#   %W% %G%
#
# Makefile for XM666
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
	../xm666/m666.f \
	../xm666/dolsq.f \
	../xm666/pm666.f 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m666

$(BINDIR)/m666: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB)  $(MDLIB) -o $(BINDIR)/m666

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m666.f
	cat $(SRCS) > ../source/m666.f

print: link
	lpr ../source/m666.f

clean:
	rm -f *.o core


