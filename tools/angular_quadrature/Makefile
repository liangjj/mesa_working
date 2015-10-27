#deck   %W% %G%
#
# Makefile for M602
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
	gauleg.f \
	lebedev.f \
	lgndrx.f \
	plgndr.f \
	radial.f \
	quad.f \
	xm602.f \
	ylm.f 
 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/xm602
$(BINDIR)/xm602: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/xm602

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/xm602.f
	cat $(SRCS) > ../source/xm602.f

print: link
	lpr ../source/xm602.f

clean:
	rm -f *.o core



