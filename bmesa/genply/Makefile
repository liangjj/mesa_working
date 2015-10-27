#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M7030
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

GET = sccs get

SRCS = \
	bdat.f \
	droppt.f \
	dropfn.f \
	genply.f \
	matcpy.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7030

$(BINDIR)/m7030: $(SRCS:.f=.o) $(PDELIB) $(DVDLIB) $(MESCTLIB) $(LAPACK) \
                               $(BLAS) $(VECLIB)  
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(PDELIB) $(DVDLIB) $(MESCTLIB) \
                                        $(LAPACK) $(BLAS) $(VECLIB) \
				-o $(BINDIR)/m7030

