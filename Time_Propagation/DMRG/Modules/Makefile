#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
#
RCS = RCS
CO = co -l
DMRG_MODULE_DIR = $(MESA_HOME)/Time_Propagation/DMRG/Modules
DMRG_MOD_LIB = $(DMRG_MODULE_DIR)/dmrg_mod_lib_f90.a
RANLIB = ranlib

MOD_SRCS = \
	$(DMRG_MODULE_DIR)/system_parameters_module.f90 \
	$(DMRG_MODULE_DIR)/Matrix_Manipulation_module.f90 \
	$(DMRG_MODULE_DIR)/Matrix_Exponential_module.f90 \
	$(DMRG_MODULE_DIR)/MPDtools_module.f90 \
	$(DMRG_MODULE_DIR)/local_operations_module.f90 \
	$(DMRG_MODULE_DIR)/io_module.f90 \
	$(DMRG_MODULE_DIR)/observables_module.f90 \
	$(DMRG_MODULE_DIR)/propagation_module.f90

.f90.o:
	$(F90) $(F90FLAGS) $<

all: $(DMRG_MOD_LIB)
$(DMRG_MOD_LIB): $(MOD_SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(DMRG_MOD_LIB) $(?:.f90=.o)
	$(RANLIB) $(DMRG_MOD_LIB)

files: $(MOD_SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod
