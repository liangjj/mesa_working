#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M206
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
	basout.f \
	m206.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m206

clean:
	rm -f *.o core *~

$(BINDIR)/m206: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m206














