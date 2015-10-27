#   @(#)Makefile	1.2   8/1/91
#
# Makefile for M390
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../m390/cblock.f \
	../m390/count.f \
	../m390/m390.f \
	../m390/matblk.f  \
	../m390/oneint.f \
	../m390/pm390.f \
	../m390/putvec.f \
	../m390/symprd.f \
	../m390/sympt.f \
	../m390/twoint.f \
	../m390/unpack.f


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m390

$(BINDIR)/m390: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/m390

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m390.f
	cat $(SRCS) > ../source/m390.f

print: link
	lpr ../source/m390.f

clean:
	rm -f *.o core


