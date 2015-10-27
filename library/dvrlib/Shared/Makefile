#   @(#)Makefile	1.4   8/2/91
#
# Makefile for dvrlib
#
DVRLIB = ../dvrlib.a
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

all: $(DVRLIB)

$(DVRLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(DVRLIB) $(?:.F=.o)

sources: $(CPPSRCS)
	co -l RCS/*.F