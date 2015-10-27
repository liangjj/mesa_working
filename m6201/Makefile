#  @(#)Makefile	1.8 8/1/91
#
# Makefile for M6201
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
REL = 

SRCS = \
	appint.f \
	bkwrd.f \
	caspln.f \
	convr.f \
	drvint.f \
	fbgint.f \
	fbnint.f \
	forwrd.f \
	ftoylm.f \
	ftoplm.f \
	legend.f \
	legint.f \
	lmdcmp.f \
	miscfn.f \
	mkspln.f \
	mkylm.f \
	nrmylm.f \
	plmtof.f \
	runiqe.f \
	scalyu.f \
	scmphi.f \
	setzro.f \
	sqrtwt.f \
	srtleg.f \
	stnzro.f \
	wave.f \
	ylmint.f \
	yuks.f  \
	yukns.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6201

$(BINDIR)/m6201: $(SRCS:.f=.o)  ${MESCTLIB} $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ${MESCTLIB} $(VECLIB) -o $(BINDIR)/m6201
