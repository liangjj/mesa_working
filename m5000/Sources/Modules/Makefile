#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/m5000
include $(TOP_DIR)/Makefile.inc
include $(MESA_HOME)/Makefile.inc
MOD_SRCS = \
	Angular_Module.f90 \
	auto_correlation_module.f90 \
	moment_module.f90 \
	plot_module.f90 \
	initial_state_module.f90 \
	Data_Subroutines_Module.f90 \
	Exponential_on_Vector_Module.f90 \
	Preconditioner_Module.f90 \
	Packed_Matrix_Module.f90 \
	Non_Packed_Matrix_Module.f90 \
	B_Spline_Matrix_Module.f90 \
	Atomic_State_Matrix_Module.f90 \
	Generalized_Lanczos_Module.f90 \
	Preconditioned_Lanczos_Module.f90 \
	Propagation_Module.f90

all: $(SIL_MOD_LIB) $(MOD_SRCS)
$(SIL_MOD_LIB): $(MOD_SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(SIL_MOD_LIB) $(?:.f90=.o)
	$(RANLIB) $(SIL_MOD_LIB)
	cp *.mod $(MODULE_DIR)/.
	cp *.o $(OBJECT_DIR)/.

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(OBJECT_DIR)/*.o
	rm -f $(MODULE_DIR)/*.mod
	rm -f $(LIB_DIR)/*.a
movein:
	mv *.mod $(MODULE_DIR)/.
	mv *.o $(OBJECT_DIR)/.

moveout:
	mv $(MODULE_DIR)/*.mod .
	mv $(OBJECT_DIR)/*.o .

