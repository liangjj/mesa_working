#   %W% %G%
#
# Makefile for XM999
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
	../xm999/m999.f \
	../xm999/dolsq.f \
	../xm999/pm999.f 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m999

$(BINDIR)/m999: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB)  $(MDLIB) -o $(BINDIR)/m999

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


