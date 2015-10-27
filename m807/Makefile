#  @(#)Makefile	5.1   11/6/94
#
# Makefile for M807
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
	../m807/m807.f \
	../m807/pm807.f \
	../m807/chkcfg.f \
	../m807/cidata.f \
	../m807/final.f \
	../m807/newcfg.f \
	../m807/btime.f \
	../m807/chorre.f \
	../m807/config.f \
	../m807/gencfg.f \
	../m807/wrtcfg.f
 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m807
$(BINDIR)/m807: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) \
		 $(MDLIB) -o $(BINDIR)/m807

sources: $(SRCS) 
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m807.f
	cat $(SRCS) > ../source/m807.f

print: link
	lpr ../source/m807.f

clean:
	rm -f *.o core


