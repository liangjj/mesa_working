#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr iterative routines
#
ITDVRLIB = ../../itdvrlib.a
BLIULIB = ../../bliulib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c -O5
GET = sccs get
REL = 

SRCS = \
	h12dvr.f \
	h1e.f \
	h12h0.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(ITDVRLIB)
$(ITDVRLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(ITDVRLIB) $(?:.f=.o)

all: $(BLIULIB)
$(BLIULIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(BLIULIB) $(?:.f=.o)
