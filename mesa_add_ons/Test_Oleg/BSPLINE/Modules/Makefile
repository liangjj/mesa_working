#   @(#)Makefile        1.4   8/2/91
#
# This Makefile will compile and make a library out of the modules needed for the BSPLINE code.
# Any main program using thedse modules needs to have the I paramete set properly in the compile line to find these
# modules and also load the library $(BSPLINE_MODULE_LIBRARY) as defined in the Makefile.inc makefile.
#
TOP_DIR = $(MESA_HOME)/Test_Oleg
include $(TOP_DIR)/Makefile.inc

SRCS = \
	spline_param.f90 \
	spline_orbitals.f90 \
	spline_atomic.f90 \
	spline_densities.f90 \
	spline_galerkin.f90 \
	spline_grid.f90 \
	spline_hl.f90 \
	spline_integrals.f90 \
	spline_moments.f90 \
	spline_Rk_integrals.f90 \
	spline_slater.f90

all: $(BSPLINE_MODULE_LIB)
$(BSPLINE_MODULE_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(BSPLINE_MODULE_LIB) $(?:.f90=.o)
	$(RANLIB) $(BSPLINE_MODULE_LIB)

clean:
	rm -f *.o
	rm -f *.mod
