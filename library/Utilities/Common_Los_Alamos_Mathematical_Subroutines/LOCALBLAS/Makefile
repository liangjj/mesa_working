#   %W%   %G%
#
#  Makefile for the LOCALLACK routines
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR = $(UTILITY)
include $(TOP_DIR)/Makefile.inc

GET = sccs get
REL = 

SRCS = \
	cgetrf.f \
        cgetrs.f \
        dgetf2.f \
        dgetrf.f \
        dgetrs.f \
        dlaswp.f \
        dlamch.f \
        dlae2.f \
        dlaev2.f \
        dlansy.f \
        dlapy2.f \
        dlartg.f \
        dlarfg.f \
        dlarft.f \
        dlarfb.f \
        dlarf.f \
        dlasr.f \
        dlassq.f \
        dlatrd.f \
        dlazro.f \
        dorgql.f \
        dorgqr.f \
        dorgtr.f \
        dorg2r.f \
        dorg2l.f \
        dsteqr.f \
        dsterf.f \
        dsyev.f \
        dsytrd.f \
        dsytd2.f \
        ilaenv.f \
        zgetf2.f \
        zgetrf.f \
        zgetrs.f \
        zlaswp.f

OBJECTS = \
	$(SRCS:.f=.o)

all: $(BLASLIB)
$(BLASLIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $(BLASLIB) $(OBJECTS)
	$(RANLIB) $(BLASLIB)
	$(AR) $(ARFLAGS) $(UTILITY_LIB) $(OBJECTS)
	$(RANLIB) $(UTILITY_LIB) 

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@
