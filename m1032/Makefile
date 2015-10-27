#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1032
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
	../m1032/cider.f \
	../m1032/m1032.f \
	../m1032/pm1032.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1032

$(BINDIR)/m1032: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1032

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1032.f
	cat $(SRCS) > ../source/m1032.f

print: link
	lpr ../source/m1032.f

clean:
	rm -f *.o core







