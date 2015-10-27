#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6298
#
FC = f77
FFLAGS = -c -O4 -recursive

LD = f77
LDFLAGS =
BINDIR = /usr1/mesa/bin
MESCTLIB = /usr1/mesa/library/mesctlib.a
GET = sccs get

SRCS = \
	filtst.f \
	recur.f \
	tstmem.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6298

$(BINDIR)/m6298: $(SRCS:.f=.o) $(MESCTLIB)  
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(MESCTLIB) -o $(BINDIR)/m6298
