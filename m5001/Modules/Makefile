#
#  Makefile for the propagator code
#
SOURCE = m5001
include $(MESA_HOME)/$(SOURCE)/Makefile.inc
SRCLOC = $(MESA_HOME)/$(SOURCE)/Modules

SRCS = \
	$(SRCLOC)/Angular_Module.f90 \
	$(SRCLOC)/auto_correlation_module.f90 \
	$(SRCLOC)/moment_module.f90 \
	$(SRCLOC)/plot_module.f90 \
	$(SRCLOC)/initial_state_module.f90 \
	$(SRCLOC)/Data_Subroutines_Module.f90 \
	$(SRCLOC)/Exponential_on_Vector_Module.f90 \
	$(SRCLOC)/Preconditioner_Module.f90 \
	$(SRCLOC)/Packed_Matrix_Module.f90 \
	$(SRCLOC)/Non_Packed_Matrix_Module.f90 \
	$(SRCLOC)/B_Spline_Matrix_Module.f90 \
	$(SRCLOC)/Atomic_State_Matrix_Module.f90 \
	$(SRCLOC)/Generalized_Lanczos_Module.f90 \
	$(SRCLOC)/Preconditioned_Lanczos_Module.f90 \
	$(SRCLOC)/Propagation_Module.f90

all: $(SIL_MOD_LIB) $(SRCS)
$(SIL_MOD_LIB): $(SRCS:.f90=.o) 
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
