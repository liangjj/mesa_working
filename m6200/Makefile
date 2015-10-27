#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6200
#
FC = f77
FFLAGS = -c -O4

LD = f77
LDFLAGS =
MAIN = /usr1/mesa
BINDIR = $(MAIN)/bin
MESCTLIB = $(MAIN)/library/mesctlib.a
LAPACK = $(MAIN)/library/LAPACK/lapack.a
BLAS = $(MAIN)/library/LAPACK/blas.a
PDELIB = $(MAIN)/library/pdelib/pdelib.a
DVDLIB = $(MAIN)/library/dvdson/dvdlib.a
VECLIB = $(MAIN)/library/libdxml.a

GET = sccs get

SRCS = \
	aij.f \
	ang.f \
	getbnd.f \
	getcnd.f \
	getdnd.f \
	gnquad.f \
	lebdev.f \
	mkgr.f \
	mkvwt.f \
	mkwt.f \
	mkyunt.f \
	radquad.f \
	rdang.f \
	rdlebdv.f \
	rdpt.f \
	satshl.f \
	scalwt.f \
	shells.f \
	sumncw.f \
	voronoi.f \
	xyz.f \
	xyzw.f \
	yukawa.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6200

$(BINDIR)/m6200: $(SRCS:.f=.o)  ${MESCTLIB} $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ${MESCTLIB} $(VECLIB) -o $(BINDIR)/m6200




