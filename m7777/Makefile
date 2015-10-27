#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M7777
#
FC = fort
FFLAGS = -c -O5

LD = fort
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB = /usr/lib/libcxml.a
GET = sccs get
REL = 

SRCS = \
	../m7777/memory.f \
	../m7777/tstmem.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7777

$(BINDIR)/m7777: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m7777

