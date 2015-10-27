#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr iterative routines
#
ITDVRLIB = ../../itdvrlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c -O5
GET = sccs get
REL = 

SRCS = \
	cblslv.f \
	cfctrs.f \
	csublu.f \
	cv2rv.f \
	rv2cv.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(ITDVRLIB)
$(ITDVRLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(ITDVRLIB) $(?:.f=.o)



