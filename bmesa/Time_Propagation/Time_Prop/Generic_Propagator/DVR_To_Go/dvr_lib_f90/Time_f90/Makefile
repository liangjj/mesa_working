#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
DVR_LIB = ../../dvr_lib_f90.a
PROP_MOD_DIR =../../../Prop_To_Go/Modules
MOD_DIR = ../../../Modules
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

CO = -co -l
RCS = RCS
REL = 

SRCS = \
	fdf.f90 \
	nrmlze.f90 \
	timpts.f90 \
	tlobato.f90 \
	tphy.f90 \
	tplmat.f90 \
	twpt.f90
.f90.o:
	$(FC) $(FFLAGS) -I$(MOD_DIR) $<

all: $(DVR_LIB)
$(DVR_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(DVR_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90*

clean:
	rm -f *.o

