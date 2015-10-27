#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the davidson eigenvalue and linear system routines
#
BLIULIB = ../bliulib.a
ITDVRLIB =../itdvrlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
CLEAN = rm

CPPFLAGS = -DDECPOINTER
FFLAGS = -c -O5 $(CPPFLAGS)
REL = 

CPPSRCS = \
	diagnl.F \
	frmres.F \
	genvec.F \
	hsmall.F \
	init.F \
	invec.F \
	outvec.F \
	ntrial.F \
	preig.F \
	rdciv.F \
	rdiag.F \
	trilin.F \
	tstovl.F

.F.o:
	$(FC) $(FFLAGS) $<

all: $(BLIULIB) $(ITDVRLIB)
$(BLIULIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(BLIULIB) $(?:.F=.o)
$(ITDVRLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(ITDVRLIB) $(?:.F=.o)

sources: $(CPPSRCS)
	co -l RCS/*.F
