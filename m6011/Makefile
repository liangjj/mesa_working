#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6011
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a
GET = sccs get

SRCS = \
	srt.f \
	srtorb.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6011

$(BINDIR)/m6011: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6011



