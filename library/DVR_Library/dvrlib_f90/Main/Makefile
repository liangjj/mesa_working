#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR=$(DVR)
include $(TOP_DIR)/Makefile.inc
# The Main Directory

RCS = RCS
CO = co -l
SRCS = \
	atomic_input.f90 \
	dvr_input.f90 \
	dvr_basis.f90 \
	ptcal.f90 \
	read_fourier.f90 \
	read_legendre.f90 \
	read_hermite.f90 \
	read_laguerre.f90 \
	read_grid_parameters.f90 \
	read_theta.f90 

all: $(DVR_LIB)
$(DVR_LIB): $(SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(DVR_LIB) $(?:.f90=.o)
	$(RANLIB) $(DVR_LIB)


sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o
