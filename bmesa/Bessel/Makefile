#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M8030
#
FC = f77
FFLAGS = -c -O4

LD = f77
LDFLAGS =
LIBDIR = /usr1/mesa/library
BINDIR = /usr1/mesa/bin
MESCTLIB = $(LIBDIR)/mesctlib.a
VECLIB = $(LIBDIR)/libdxml.a

GET = sccs get

SRCS = \
	drive.f \

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m8030

$(BINDIR)/m8030: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)  
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m8030

