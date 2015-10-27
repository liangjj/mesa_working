#   @(#)Makefile	1.9   8/7/91
#
#  
#
TOP_DIR = $(MESA_LIB)/Mesalib_f90
include $(TOP_DIR)/Makefile.inc
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
	FEDVR_Shared.f90 \
	FEDVR_Derived_Types.f90 \
	Matrix_Print.f90 

all: $(MODULE_LIB)  
$(MODULE_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(MODULE_LIB) $(?:.f90=.o)
	$(AR) $(ARFLAGS) $(UTILITY_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod

