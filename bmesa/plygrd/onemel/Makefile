#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6291
#
FC = f77
FFLAGS = -c -O4

LD = f77
LDFLAGS =
BINDIR = /usr1/mesa/bin
MESCTLIB = /usr1/mesa/library/mesctlib.a
LAPACK = /usr1/mesa/library//LAPACK/lapack.a
BLAS = /usr1/mesa/library/LAPACK/blas.a
PDELIB = /usr1/mesa/library/pdelib/pdelib.a
DVDLIB = /usr1/mesa/library/dvdson/dvdlib.a
VECLIB = /usr1/mesa/library/libdxml.a

#BINDIR = /usr1/users/bis/mesa/bin
#MESCTLIB = /usr1/users/bis/mesa/library/mesctlib.a
#LAPACK = /usr1/users/bis/mesa/library//LAPACK/lapack.a
#BLAS = /usr1/users/bis/mesa/library/LAPACK/blas.a
#PDELIB = /usr1/users/bis/mesa/library/pdelib/pdelib.a
#DVDLIB = /usr1/users/bis/mesa/library/dvdson/dvdlib.a
#VECLIB = /usr1/users/bis/mesa/library/libdxml.a

GET = sccs get

SRCS = \
	addd.f \
	h0.f \
	h1.f \
	one.f \
	tocmplx.f \
	vmat1d.f \
	vt0.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6291

$(BINDIR)/m6291: $(SRCS:.f=.o) $(PDELIB) $(DVDLIB) $(MESCTLIB) $(LAPACK) \
                               $(BLAS) $(VECLIB)  
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(PDELIB) $(DVDLIB) $(MESCTLIB) \
                                        $(LAPACK) $(BLAS) $(VECLIB) \
					-o $(BINDIR)/m6291
