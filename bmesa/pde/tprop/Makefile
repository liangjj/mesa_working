#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M7080
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
TDIFFLIB = /usr1/mesa/bmesa/pde/tprop/tdifflib/tdifflib.a

GET = sccs get

SRCS = \
	cfactr.f \
	derivs.f \
	df.f \
	gpaket.f \
	grid.f \
	h0.f \
	honv.f \
	prdiff.f \
	rfactr.f \
	setind.f \
	srteig.f \
	tdrver.f \
	tprop.f \
	vpert.f \
	wrdadd.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7080

$(BINDIR)/m7080: $(SRCS:.f=.o) $(TDIFFLIB) $(PDELIB) $(DVDLIB) $(MESCTLIB) \
                               $(LAPACK)  $(BLAS) $(LIBDXML)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(TDIFFLIB) $(PDELIB) $(DVDLIB) \
                                       $(MESCTLIB) $(LAPACK) $(BLAS) \
                                       $(LIBDXML) -o $(BINDIR)/m7080
