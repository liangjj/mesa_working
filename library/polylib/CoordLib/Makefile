#   @(#)Makefile	1.4   8/2/91
#
# Makefile for coordlib
#
COORDLIB = ../coordlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
CPPFLAGS = -DDECPOINTER
FFLAGS = -c -O4 $(CPPFLAGS)

CPPSRCS = \
	cpoly.F \
	drvply.F \
	lgngr.F \
	cordfn.F

.F.o:
	$(FC) $(FFLAGS) $<

all: $(COORDLIB)

$(COORDLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(COORDLIB) $(?:.F=.o)

sources: $(CPPSRCS)
	co -l RCS/*.F