#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M7101
#
FC = fort
FFLAGS = -c -O4

LD = fort
LDFLAGS =

MAIN = /usr/home/bis/mesa
BINDIR = $(MAIN)/bin
DVRLIB = $(MAIN)/library/dvrlib/dvrlib.a
POTLIB = $(MAIN)/library/potlib/potlib.a
HAMLIB = $(MAIN)/bmesa/R_Matrix/hamlib/hamlib.a
MESALIB = $(MAIN)/library/mesalib.a
PDELIB = $(MAIN)/library/pdelib/pdelib.a
DVDLIB = $(MAIN)/library/dvdson/dvdlib.a
MDLIB = /usr/lib/libcxml.a


GET = sccs get

SRCS = \
	bnd1d.f \
	conrmat.f \
	dvddat.f \
	fobj.f \
	gam1.f \
	gtrbox.f \
	gkqunsy.f \
	gkqsy.f \
	gusums.f \
	gssums.f \
	hamio.f \
	h2h0.f \
	intrv.f \
	lobato.f \
	lindat.f \
	mke.f \
	nwtrap.f \
	ovsym.f \
	ovusym.f \
	ov1d.f \
	rmtrx1.f \
	rgamma.f \
	vmat.f


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7101

$(BINDIR)/m7101: $(SRCS:.f=.o) $(DVRLIB) $(HAMLIB) $(POTLIB) $(DVDLIB) $(PDELIB) \
                               $(MESALIB) $(MDLIB)    
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(DVRLIB) $(HAMLIB) $(POTLIB) $(DVDLIB) \
                                        $(MESALIB) $(PDELIB) $(MDLIB) \
					-o $(BINDIR)/m7101

