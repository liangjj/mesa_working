#  @(#)Makefile	1.8 8/1/91
#
# Makefile for M6203
#
FC = f77
FFLAGS = -c  

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
	ftoplm.f \
	ftoylm.f \
	homog.f \
	legend.f \
	legint.f \
	lmdcmp.f \
	miscfn.f \
	mkspln.f \
	mkylm.f \
	nrmylm.f \
	plmtof.f \
	poisson.f \
	runiqe.f \
	scalyu.f \
	scmphi.f \
	setzro.f \
	sqrtwt.f \
	srtleg.f \
	ylmint.f \
	yukns.f \
	yuks.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6203

$(BINDIR)/m6203: $(SRCS:.f=.o)  ${MESCTLIB} $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ${MESCTLIB} $(VECLIB) -o $(BINDIR)/m6203
