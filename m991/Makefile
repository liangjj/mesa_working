

#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M991
#
FC = fort 
FFLAGS = -c 
LD = fort
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB = /usr/lib/libcxml.a
GET = sccs get
REL = 

SRCS = \
	m991.f \
	pm991.f \
	iomesa.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m991

$(BINDIR)/m991: $(SRCS:.f=.o) $(MESALIB) $(MDLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m991

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m991.f
	cat $(SRCS) > ../source/m991.f

print: link
	lpr ../source/m991.f

clean:
	rm -f *.o core



