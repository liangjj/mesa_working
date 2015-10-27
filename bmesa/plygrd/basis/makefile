#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6290
#
FC = f77
FFLAGS = -c -O4

LD = f77
LDFLAGS =
LIBDIR = /usr1/mesa/library
BINDIR = /usr1/mesa/bin
MESCTLIB = $(LIBDIR)/mesctlib.a
PDELIB = $(LIBDIR)/pdelib/pdelib.a
VECLIB = $(LIBDIR)/libdxml.a

GET = sccs get

SRCS = \
	basis.f \
	boxfn.f \
	cheby.f \
	chebpt.f \
	cmpare.f \
	dropfn.f \
	droppt.f \
	cfine.f \
	expand.f \
	gdat.f \
	gi2gj.f \
	gwadd.f \
	i2j1d.f \
	i2j2d.f \
	i2j3d.f \
	i2j4d.f \
	prepfn.f \
	tdat.f \
	v1di2o.f \
	v2di2o.f \
	v3di2o.f \
	v4di2o.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6290

$(BINDIR)/m6290: $(SRCS:.f=.o) $(PDELIB) $(MESCTLIB)   
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(PDELIB) $(MESCTLIB) \
					-o $(BINDIR)/m6290

