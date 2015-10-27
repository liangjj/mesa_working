#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Atomic_Propagator
include $(TOP_DIR)/Makefile.inc
MOD_SRCS = \
	Angular_Module.f90 \
	auto_correlation_module.f90 \
	moment_module.f90 \
	plot_module.f90 \
	initial_state_module.f90 \
	spatial_wavefunction_module.f90 \
	Data_Subroutines_Module.f90 \
	Preconditioner_Module.f90 \
	Packed_Matrix_Module.f90 \
	B_Spline_Matrix_Module.f90 \
	Atomic_State_Matrix_Module.f90 \
	Exponential_on_Vector_Module.f90 \
	Lanczos_Module.f90 \
	Preconditioned_Lanczos_Module.f90 \
	Propagation_Module.f90



.mod.o:
	
#.f90.o:
#	$(F90) $(F90FLAGS) -I$(MODULES) $<
#.f90.mod:
#	$(F90) $(F90FLAGS) -I$(MODULES) $<

all: $(ATOMIC_PROP_MOD_LIB)
$(ATOMIC_PROP_MOD_LIB): $(MOD_SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(ATOMIC_PROP_MOD_LIB) $(?:.f90=.o)
	$(RANLIB) $(ATOMIC_PROP_MOD_LIB)

files: $(MOD_SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
