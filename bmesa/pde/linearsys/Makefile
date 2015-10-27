#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M7040
#
FC = f77
FFLAGS = -c +OP4  

LD = f77
LDFLAGS =
BINDIR = /newdsk/mesa/bin
MESCTLIB = /newdsk/mesa/library/mesctlib.a
LAPACK = /newdsk/mesa/library//LAPACK/lapack.a
BLAS = /newdsk/mesa/library/LAPACK/blas.a
DVDLIB = /newdsk/mesa/library/dvdson/dvdlib.a
PDELIB = /newdsk/mesa/library/pdelib/pdelib.a
VECLIB = /usr/lib/libvec.a

GET = sccs get

SRCS = \
	calrhs.f \
	exact.f \
	genmat.f \
	lindat.f \
	lschr.f \
	pde.f \
	tocord.f 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7040

$(BINDIR)/m7040: $(SRCS:.f=.o) $(PDELIB) $(DVDLIB) $(MESCTLIB) $(LAPACK) \
                               $(BLAS) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(PDELIB) $(DVDLIB) $(MESCTLIB) $(LAPACK) \
                                       $(BLAS) $(VECLIB) \
	                                -o $(BINDIR)/m7040
