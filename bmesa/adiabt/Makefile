#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M7100
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
	adiabt.f \
	fileig.f \
	rdlabl.f \
	setham.f \
	vparm.f \
	vadiab.f \
	wrtfil.f \
	ylamda.f


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7100

$(BINDIR)/m7100: $(SRCS:.f=.o) $(DVRLIB) $(HAMLIB) $(POTLIB) $(MESCTLIB) \
			       $(LAPACK) $(VECLIB)    
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(DVRLIB) $(HAMLIB) $(POTLIB) \
					$(MESCTLIB) $(LAPACK) $(VECLIB) \
					-o $(BINDIR)/m7100
clean:	
	rm -f *.o
