#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M8005
#
FC = f77
FFLAGS = -c -O4

LD = f77
LDFLAGS =
LIBDIR = /usr1/mesa/library
BINDIR = /usr1/mesa/bin
DVRLIB = $(LIBDIR)/dvrlib/dvrlib.a
TDVRLIB = $(LIBDIR)/tdvrlib/tdvrlib.a
POTLIB = $(LIBDIR)/potlib/potlib.a
MESCTLIB = $(LIBDIR)/mesctlib.a
LAPACK = $(LIBDIR)/LAPACK/lapack.a
BLAS = $(LIBDIR)/LAPACK/blas.a
PDELIB = $(LIBDIR)/pdelib/pdelib.a
DVDLIB = $(LIBDIR)/dvdson/dvdlib.a
ITDVRLIB = $(LIBDIR)/itdvr/itdvrlib.a
VECLIB = $(LIBDIR)/libdxml.a

GET = sccs get

SRCS = \
	cctdvr.f \
	psit.f 


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m8005

$(BINDIR)/m8005: $(SRCS:.f=.o) $(TDVRLIB) $(DVRLIB) $(POTLIB) $(ITDVRLIB) \
                               $(PDELIB) $(DVDLIB)  $(MESCTLIB) $(LAPACK) \
                               $(BLAS) $(VECLIB)  
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(TDVRLIB) $(DVRLIB) $(POTLIB) $(ITDVRLIB) \
                                        $(PDELIB) $(DVDLIB) $(MESCTLIB) $(LAPACK) \
                                        $(BLAS) $(VECLIB) -o $(BINDIR)/m8005

