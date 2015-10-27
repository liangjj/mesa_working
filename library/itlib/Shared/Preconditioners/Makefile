#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the preconditioners library
#
BLIULIB = ../../bliulib.a
ITDVRLIB =../../itdvrlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c -O5
GET = sccs get
REL = 

SRCS = \
	blkout.f \
	blkpre.f \
	diddle.f \
	dslve.f \
	dcslve.f \
	escale.f \
	h0prec.f \
	h0cprec.f \
	rbldag.f \
	rblslv.f \
	rblvec.f \
	rfctrs.f \
	rsublu.f \
	seppre.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BLIULIB) $(ITDVRLIB)
$(BLIULIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(BLIULIB) $(?:.f=.o)
$(ITDVRLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(ITDVRLIB) $(?:.f=.o)





