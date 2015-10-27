#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)

AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l
MOD_SRCS = \
	auto_correlation_module.f90 \
	moment_module.f90 \
	plot_module.f90 \
	regional_module.f90 \
	normalize_module.f90 \
	non_linear_potential_module.f90 \
	initial_state_module.f90 \
	check_convergence_module.f90 \
	check_normalization_module.f90 \
	test_overlap_module.f90 \
	spatial_wavefunction_module.f90 \
	matrix_vector_multiply_module.f90 \
	arnoldi_module.f90 \
	split_operator_module.f90 

%.mod : %.o
	@if [! -f $@ ]; then \
	rm $< \
	$(MAKE) $< \
	fi

%.o : %.f90
	$(F90) $(F90FLAGS) -I$(MOD_DIR) $<

all: $(MESA_PROPAGATOR_LIB)
$(MESA_PROPAGATOR_LIB): $(MOD_SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(MESA_PROPAGATOR_LIB) $(?:.f90=.o)
	$(RANLIB) $(MESA_PROPAGATOR_LIB)

files: $(MOD_SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o 
	rm -f *.mod
