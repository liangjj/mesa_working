#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6290
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
	ainver.f \
	bdat.f \
	boxfn.f \
	calfun.f \
	cheby.f \
	chebpt.f \
	cmpare.f \
	cmpre.f \
	cmpslv.f \
	cor2cor.f \
	dropfn.f \
	droppt.f \
	ccorse.f \
	cfine.f \
	expand.f \
	fildr1.f \
	fildr2.f \
	fildr3.f \
	fin2cor.f \
	gdat.f \
	gpaket.f \
	gwadd.f \
	hamfl1.f \
	hamfl2.f \
	hamfl3.f \
	hamtot.f \
	matcpy.f \
	mkh1d.f \
	mkh2d.f \
	mkh3d.f \
	mkpsi0.f \
	mkhpsi.f \
	prepfn.f \
	setd2.f \
	setd3.f \
	setd4.f \
	setnd2.f \
	setnd3.f \
	setnd4.f \
	vect0.f \
	wrd.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6290

$(BINDIR)/m6290: $(SRCS:.f=.o) $(PDELIB) $(DVDLIB) $(MESCTLIB) $(LAPACK) \
                               $(BLAS) $(VECLIB)  
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(PDELIB) $(DVDLIB) $(MESCTLIB) \
                                        $(LAPACK) $(BLAS) $(VECLIB) \
					-o $(BINDIR)/m6290

