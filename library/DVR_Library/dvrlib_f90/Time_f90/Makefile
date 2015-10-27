#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
.SUFFIXES: .f90 $(SUFFIXES)
DVR_LIB = ../../dvr_lib_f90.a
TOP_DIR=$(DVR)
include $(TOP_DIR)/Makefile.inc
RCS = RCS
CO = co -l
REL = 
RANLIB = ranlib

SRCS = \
	fdf.f90 \
	nrmlze.f90 \
	timpts.f90 \
	tlobato.f90 \
	tphy.f90 \
	tplmat.f90 \
	twpt.f90

all: $(DVR_LIB)
$(DVR_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(DVR_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90*

clean:
	rm -f *.o

