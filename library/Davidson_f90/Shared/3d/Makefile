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
	h32dvr.f \
	h3e.f \
	h32h0.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(ITDVRLIB)
$(ITDVRLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(ITDVRLIB) $(?:.f=.o)

all: $(BLIULIB)
$(BLIULIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(BLIULIB) $(?:.f=.o)
