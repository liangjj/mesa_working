#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
include $(MODULES)/Makefile.inc
.SUFFIXES: .f90 .mod $(SUFFIXES)

SRCS = \
	accuracy.f90 \
	input_output.f90 \
	Data_Module.f90 \
	potential.f90 \
	grid_global.f90 \
	dvr_shared.f90 \
	dvr_global.f90 \
	Atomic_Matrices.f90 \
	Global_Time_Propagation_Module.f90 \
	prop_prnt.f90 \
	prop_global.f90 \
	arnoldi_global.f90\
	lanczos_global.f90\
	lanczos_prnt.f90\
	arnoldi_global_it.f90 \
	arnoldi_global_rt.f90 \
	dvd_prnt.f90 \
	dvr_prnt.f90 \
	dvr_global.f90 \
	fd_prnt.f90 \
	fd_global.f90 \
	dvrprop_global_rt.f90 \
	dvrprop_global_it.f90 \
	dvrprop_global.f90 \
	Pack_Global.f90 \
	Iterative_Global.f90 \
	Pack_Hamiltonian_Module.f90 \
	Matrix_Utility_Module.f90 \
	r_matrix_global.f90 \
	mesa_global.f90 \
	dvr_matrix_vector_multiply_module.f90 \
	full_matrix_vector_multiply_module.f90 \
	full_matrix_vector_iteration_module.f90 \
	packed_matrix_vector_multiply_module.f90 \
	Legendre.f90 \
	Legendre_Data.f90 \
	FEDVR_Global.f90 \
	FEDVR_Shared.f90 \
	Two_Electron_Shared.f90 \
	Coulomb_Shared.f90 \
	FEDVR_Derived_Types.f90 \
	Matrix_Defined_Types.f90
OBJECTS = \
	$(SRCS:.f90=.o) 
#
all: $(GEN_MOD_LIB)  
$(GEN_MOD_LIB): $(OBJECTS)  
	$(AR) $(ARFLAGS) $(GEN_MOD_LIB) $(OBJECTS)
	$(RANLIB) $(GEN_MOD_LIB)
	$(AR) $(ARFLAGS) $(MOD_LIB) $(OBJECTS)
	$(RANLIB) $(MOD_LIB)
sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod
