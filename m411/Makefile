#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M411
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
	../m411/aeqbc.f \
	../m411/cpfrez.f \
	../m411/linvec.f \
	../m411/m411.f \
	../m411/pm411.f \
	../m411/putsym.f \
	../m411/schmab.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m411

$(BINDIR)/m411: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m411

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m411.f
	cat $(SRCS) > ../source/m411.f

print: link
	lpr ../source/m411.f

clean:
	rm -f *.o core


