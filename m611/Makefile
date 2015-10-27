#   @(#)Makefile	5.2 2/5/95
#
# Makefile for M611
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
	../m611/charges.f \
	../m611/cubic.f \
	../m611/m611.f \
	../m611/pm611.f \
	../m611/sizes.f \
	../m611/spline3.f \
	../m611/splint.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m611

$(BINDIR)/m611: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB)  $(MDLIB) -o $(BINDIR)/m611

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m611.f
	cat $(SRCS) > ../source/m611.f

print: link
	lpr ../source/m611.f

clean:
	rm -f *.o core


