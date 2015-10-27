#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M904
#
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

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
	../m904/m904.f \
	../m904/pm904.f \
	../m904/analyz.f \
	../m904/btime.f \
	../m904/ciout.f \
	../m904/detdrv.f \
	../m904/detin.f \
	../m904/detout.f \
	../m904/diagnl.f \
	../m904/geneig.f \
	../m904/hinout.f \
	../m904/hrdwrt.f \
	../m904/indexf.f \
	../m904/msedet.f \
	../m904/msfdet.f \
        ../m904/msidet.f \
	../m904/msodet.f \
	../m904/oneout.f \
	../m904/output.f \
	../m904/pck.f \
	../m904/rdints.f \
	../m904/sefcig.f \
	../m904/simeig.f \
	../m904/socig.f \
	../m904/sort.f \
	../m904/spinci.f \
	../m904/sumphz.f \
	../m904/sums.f \
	../m904/twoout.f \
	../m904/soanal.f \
	../m904/unpck.f
 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m904
$(BINDIR)/m904: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) \
		 $(MDLIB) -o $(BINDIR)/m904

sources: $(SRCS) 
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m904.f
	cat $(SRCS) > ../source/m904.f

print: link
	lpr ../source/m904.f

clean:
	rm -f *.o core


