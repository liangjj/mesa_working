#   @(#)Makefile	5.1   11/6/94
#
# Makefile for tst
#
FC = fort
FFLAGS = -c

LD = fort
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB = /usr/lib/libcxml.a
REL = 
GET =

SRCS = \
	../tst/tst.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/tst

$(BINDIR)/tst: $(SRCS:.f=.o) $(MESALIB) $(MDLIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/tst

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

