#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Lanczos_Library
include $(TOP_DIR)/Makefile.inc
MOD_SRCS = \
	initial_state_module.f90 \
	Pack_Matrix_Module.f90 \
	Preconditioner_Module.f90 \
	Matrix_Module.f90 \
	Lanczos_Module.f90 \
	Lanczos_Driver_Module.f90

.mod.o:
	
.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) $<
.f90.mod:
	$(F90) $(F90FLAGS) -I$(MODULES) $<

all: $(LANCZOS_MOD_LIB)
$(LANCZOS_MOD_LIB): $(MOD_SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(LANCZOS_MOD_LIB) $(?:.f90=.o)
	$(RANLIB) $(LANCZOS_MOD_LIB)

files: $(MOD_SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
