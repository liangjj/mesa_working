#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M812
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
	../m812/check.f \
	../m812/chknrm.f \
	../m812/m812.f \
	../m812/pm812.f \
	../m812/trn1dm.f \
	../m812/trn2dm.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m812

$(BINDIR)/m812: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m812

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m812.f
	cat $(SRCS) > ../source/m812.f

print: link
	lpr ../source/m812.f

clean:
	rm -f *.o core



