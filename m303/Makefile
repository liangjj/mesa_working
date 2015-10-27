#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M303
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
	../m303/fmoned.f \
	../m303/fmt.f \
	../m303/m303.f \
	../m303/oneint.f \
	../m303/pm303.f \
	../m303/sints.f \
	../m303/tints.f \
	../m303/tranhs.f \
	../m303/transd.f \
	../m303/vints.f \
	../m303/lpints.f \
	../m303/lpd1.f \
	../m303/ltab.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m303

$(BINDIR)/m303: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m303

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m303.f
	cat $(SRCS) > ../source/m303.f

print: link
	lpr ../source/m303.f

clean:
	rm -f *.o core



