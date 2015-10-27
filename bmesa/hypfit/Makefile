#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6299
#
FC = f77
FFLAGS = -c -O4

LD = f77
LDFLAGS =

MAIN = /usr1/mesa
BINDIR = $(MAIN)/bin
DVRLIB = $(MAIN)/library/dvrlib/dvrlib.a
POTLIB = $(MAIN)/library/potlib/potlib.a
HAMLIB = $(MAIN)/library/hamlib/hamlib.a
MESCTLIB = $(MAIN)/library/mesctlib.a
LAPACK = $(MAIN)/library/LAPACK/lapack.a
BLAS = $(MAIN)/library/LAPACK/blas.a
VECLIB = $(MAIN)/library/libdxml.a


GET = sccs get

SRCS = \
	apdr1r2.f \
	cartder.f \
	cmphy.f \
	cmpr1r2.f \
	coefs.f \
	derang.f \
	derr.f \
	drvhyp.f \
	exdr1r2.f \
	fileig.f \
	hypder.f \
	hypang.f \
	hypfit.f \
	hyprad.f \
	mang.f \
	prodfn.f \
	prnfn.f \
	prntdr.f \
	radang.f \
	radial.f


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6299

$(BINDIR)/m6299: $(SRCS:.f=.o) $(DVRLIB) $(HAMLIB) $(POTLIB) $(MESCTLIB) \
                               $(LAPACK) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(DVRLIB) $(HAMLIB) $(POTLIB) \
                                        $(MESCTLIB) $(LAPACK) $(VECLIB) \
					-o $(BINDIR)/m6299
clean:	
	rm -f *.o
