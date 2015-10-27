#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1031
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
	../m1031/cilag.f \
	../m1031/m1031.f \
	../m1031/pm1031.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1031

$(BINDIR)/m1031: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1031

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1031.f
	cat $(SRCS) > ../source/m1031.f

print: link
	lpr ../source/m1031.f

clean:
	rm -f *.o core







