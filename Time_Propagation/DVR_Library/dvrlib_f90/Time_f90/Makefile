#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
.SUFFIXES: .f90 $(SUFFIXES)
MOD_DIR = $(MODULES)
DVR_LIB = $(DVR_LIBRARY_G)/dvr_lib_f90.a
RCS = RCS
CO = co -l
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
	$(F90) $(F90FLAGS) -I$(MOD_DIR) $<

all: $(DVR_LIB)
$(DVR_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(DVR_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90*

clean:
	rm -f *.o

