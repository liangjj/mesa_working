#   @(#)Makefile	1.3   8/1/91
#
# Makefile for lstsq 
#
FC = f77
FFLAGS = -c 

LD = f77
LDFLAGS =

BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a
GET = sccs get
REL = 

SRCS = \
	../m6026/lstsq.f \
	../m6026/tstlsq.f


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6026

$(BINDIR)/m6026: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6026

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m6026.f
	cat $(SRCS) > ../source/m6026.f

print: link
	lpr ../source/m101.f

clean:
	rm -f *.o core
