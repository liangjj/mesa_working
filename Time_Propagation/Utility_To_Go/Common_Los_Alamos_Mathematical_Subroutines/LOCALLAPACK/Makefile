#   %W%   %G%
#
#  Makefile for the LOCALLACK routines
#
LAPACKLIB = ../../lapacklib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c
GET = sccs get
REL = 

SRCS = \
	../LOCALLAPACK/cgetrf.f \
        ../LOCALLAPACK/cgetrs.f \
        ../LOCALLAPACK/dgetf2.f \
        ../LOCALLAPACK/dgetrf.f \
        ../LOCALLAPACK/dgetrs.f \
        ../LOCALLAPACK/dlaswp.f \
        ../LOCALLAPACK/dlamch.f \
        ../LOCALLAPACK/dlae2.f \
        ../LOCALLAPACK/dlaev2.f \
        ../LOCALLAPACK/dlansy.f \
        ../LOCALLAPACK/dlapy2.f \
        ../LOCALLAPACK/dlartg.f \
        ../LOCALLAPACK/dlarfg.f \
        ../LOCALLAPACK/dlarft.f \
        ../LOCALLAPACK/dlarfb.f \
        ../LOCALLAPACK/dlarf.f \
        ../LOCALLAPACK/dlasr.f \
        ../LOCALLAPACK/dlassq.f \
        ../LOCALLAPACK/dlatrd.f \
        ../LOCALLAPACK/dlazro.f \
        ../LOCALLAPACK/dorgql.f \
        ../LOCALLAPACK/dorgqr.f \
        ../LOCALLAPACK/dorgtr.f \
        ../LOCALLAPACK/dorg2r.f \
        ../LOCALLAPACK/dorg2l.f \
        ../LOCALLAPACK/dsteqr.f \
        ../LOCALLAPACK/dsterf.f \
        ../LOCALLAPACK/dsyev.f \
        ../LOCALLAPACK/dsytrd.f \
        ../LOCALLAPACK/dsytd2.f \
        ../LOCALLAPACK/ilaenv.f \
        ../LOCALLAPACK/zgetf2.f \
        ../LOCALLAPACK/zgetrf.f \
        ../LOCALLAPACK/zgetrs.f \
        ../LOCALLAPACK/zlaswp.f

.f.o:
	$(F90) $(F90FLAGS) $<

all: $(LAPACKLIB)

$(LAPACKLIB): $(SRCS:.f=.o) $(FOBJ)
	$(AR) $(ARFLAGS) $(LAPACKLIB) $(?:.f=.o)
	$(RANLIB) $(LAPACKLIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@
