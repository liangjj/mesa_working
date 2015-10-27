#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M7020
#
FC = f77
FFLAGS = -c   

LD = f77
LDFLAGS =
BINDIR = /usr1/mesa/bin
MESCTLIB = /usr1/mesa/library/mesctlib.a
LAPACK = /usr1/mesa/library//LAPACK/lapack.a
BLAS = /usr1/mesa/library/LAPACK/blas.a
DVDLIB = /usr1/mesa/library/dvdson/dvdlib.a
PDELIB = /usr1/mesa/library/pdelib/pdelib.a
LIBDXML =/usr1/mesa/library/libdxml.a
VECLIB =/usr1/mesa/library/mdutil/dec/veclib/veclib.a

GET = sccs get

SRCS = \
	dvddat.f \
	genmat.f \
	indset.f \
	lschr.f \
	pde.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7020

$(BINDIR)/m7020: $(SRCS:.f=.o) $(PDELIB) $(DVDLIB) $(MESCTLIB) $(LAPACK) \
                               $(BLAS) $(VECLIB) $(LIBDXML)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(PDELIB) $(DVDLIB) $(MESCTLIB) \
                                      $(LAPACK) $(BLAS) $(VECLIB) $(LIBDXML) \
                                     -o $(BINDIR)/m7020
