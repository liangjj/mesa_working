#   @(#)Makefile	1.4   8/2/91
#
# Makefile for orthlib
#
ORTHLIB = ../orthlib.a
COORDLIB = ../coordlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
CPPFLAGS = -DDECPOINTER
FFLAGS = -c -O4 $(CPPFLAGS)

CPPSRCS = \
	blcyl.F \
	blrad.F \
	blxyz.F \
	cnvtab.F \
	cnvtpt.F \
	genwtf.F \
	genq.F \
	genrwt.F \
	kemel.F \
	kinxyz.F \
	kinrad.F \
	kincyl.F \
	lancz.F \
	modab.F \
	plywpt.F \
	recur.F \
	scaprd.F \
	todelta.F \
	xsq2x.F 

.F.o:
	$(FC) $(FFLAGS) $<

all: $(ORTHLIB) $(COORDLIB)

$(ORTHLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(ORTHLIB) $(?:.F=.o)

$(COORDLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(COORDLIB) $(?:.F=.o)

sources: $(CPPSRCS)
	co -l RCS/*.F