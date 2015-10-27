#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M7103
#
FC = f77
FFLAGS = -c -O4

LD = f77
LDFLAGS =

MAIN = /usr1/mesa
#MAIN = /usr1/users/bis/mesa
BINDIR = $(MAIN)/bin
DVRLIB = $(MAIN)/library/dvrlib/dvrlib.a
POTLIB = $(MAIN)/library/potlib/potlib.a
HAMLIB = $(MAIN)/library/hamlib/hamlib.a
MESCTLIB = $(MAIN)/library/mesctlib.a
LAPACK = $(MAIN)/library/LAPACK/lapack.a
BLAS = $(MAIN)/library/LAPACK/blas.a
PDELIB = $(MAIN)/library/pdelib/pdelib.a
DVDLIB = $(MAIN)/library/dvdson/dvdlib.a
VECLIB = $(MAIN)/library/libdxml.a


GET = sccs get

SRCS = \
	drham.f \
	dvddat.f \
	ham2d.f \
	hypsph.f \
	hypfn.f \
	ind2d.f \
	itham.f \
	lindat.f \
	rdlabl.f \
	setham.f \
	vadiab.f \
	vparm.f \
	vcouhy.f \
	vexphy.f \
	vhmohy.f \
	vwelhy.f 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7103

$(BINDIR)/m7103: $(SRCS:.f=.o) $(DVRLIB) $(HAMLIB) $(POTLIB) $(DVDLIB) $(PDELIB) \
                               $(MESCTLIB) $(LAPACK) $(VECLIB)    
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(DVRLIB) $(HAMLIB) $(POTLIB) $(DVDLIB) \
                                        $(MESCTLIB) $(PDELIB) $(LAPACK) $(VECLIB) \
					-o $(BINDIR)/m7103

