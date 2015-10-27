#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
.SUFFIXES: .f90 $(SUFFIXES)
GENERIC_PROPAGATOR_DIR = \
$(MESA_HOME)/bmesa/Time_Propagation/Time_Prop/Generic_Propagator
DVR_LIBRARY =  \
$(GENERIC_PROPAGATOR_DIR)/DVR_Library
DVR_LIB_DIR = $(DVR_LIBRARY)/dvrlib_f90
FD_LIB_DIR = $(DVR_LIBRARY)/fdlib_f90
MOD_DIR = $(GENERIC_PROPAGATOR_DIR)/Modules
DVR_LIB = $(DVR_LIBRARY)/dvr_lib_f90.a
FD_LIB = $(DVR_LIBRARY)/fd_lib_f90.a

AR = ar
ARFLAGS = rvu

FFLAGS = -c -O4 -I$(MOD_DIR)
LD =$(FC)
LDFLAGS =
RCS = RCS
CO = co -l
REL = 

SRCS = \
	blcyl.f90 \
	blrad.f90 \
	blxyz.f90 \
	cnvtab.f90 \
	cnvtpt.f90 \
	genwtf.f90 \
	genq.f90 \
	genrwt.f90 \
	kemel.f90 \
	kinxyz.f90 \
	kinrad.f90 \
	kincyl.f90 \
	lancz.f90 \
	modab.f90 \
	plywpt.f90 \
	recur.f90 \
	scaprd.f90 \
	todelta.f90 \
	xsq2x.f90

.f90.o:
	$(FC) $(FFLAGS) $<

all: $(DVR_LIB)
$(DVR_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(DVR_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90*

clean:
	rm *.o

