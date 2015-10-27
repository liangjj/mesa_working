#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Generic_Propagator
include $(TOP_DIR)/Makefile.inc

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

.mod.o:

.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) $<

.f90.mod:
	$(F90) $(F90FLAGS) -I$(MODULES) $<


all: $(PROP_MOD_LIB)
$(PROP_MOD_LIB): $(MOD_SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(PROP_MOD_LIB) $(?:.f90=.o)
	$(RANLIB) $(PROP_MOD_LIB)

files: $(MOD_SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.a
