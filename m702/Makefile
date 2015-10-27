#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M702
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
	../m702/fmoned.f \
	../m702/fmt.f \
	../m702/lpd1.f \
	../m702/lpints.f \
	../m702/ltab.f \
	../m702/m702.f \
	../m702/movd2e.f \
	../m702/nucrep.f \
	../m702/nuctwo.f \
	../m702/oneint.f \
	../m702/pm702.f \
	../m702/sints.f \
	../m702/tints.f \
	../m702/transd.f \
	../m702/transf.f \
	../m702/vints.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m702

$(BINDIR)/m702: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m702

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m702.f
	cat $(SRCS) > ../source/m702.f

print: link
	lpr ../source/m702.f

clean:
	rm -f *.o core


