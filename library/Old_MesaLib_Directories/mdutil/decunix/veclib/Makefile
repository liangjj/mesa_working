#   @(#)Makefile	1.7   8/7/91
#
#  Makefile for the machine-dependent portion of MesaLib
#
VECLIB = veclib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib 

FC = f77
FFLAGS = -c
GET = sccs get
REL = 

SRCS = \
	../veclib/cdotc.f \
	../veclib/cdotu.f \
	../veclib/cgemm.f \
	../veclib/saxpy.f \
	../veclib/scatter.f \
	../veclib/scopy.f \
	../veclib/sdot.f \
	../veclib/sgemm.f \
	../veclib/sgemv.f \
	../veclib/sscal.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(VECLIB) 
$(VECLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(VECLIB) $(?:.f=.o)
	$(RANLIB) $(VECLIB)
