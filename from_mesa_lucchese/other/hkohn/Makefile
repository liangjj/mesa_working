#   %W%   %G%
#
# Makefile for hkohn
#

FC = g77
FFLAGS = -finit-local-zero -c 

LD = g77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
GET = sccs get
REL = 

SRCS = \
	./hkohn.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/xhkohn

$(BINDIR)/xhkohn: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/xhkohn




