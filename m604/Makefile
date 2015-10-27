#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M604
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
GET = sccs get
REL = 

SRCS = \
	../m604/detsym.f \
	../m604/m604.f \
	../m604/matprt.f \
	../m604/nsalc.f \
	../m604/pm604.f \
	../m604/symorb.f \
	../m604/tstsym.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m604

$(BINDIR)/m604: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m604

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m604.f
	cat $(SRCS) > ../source/m604.f

print: link
	lpr ../source/m604.f

clean:
	rm -f *.o core


