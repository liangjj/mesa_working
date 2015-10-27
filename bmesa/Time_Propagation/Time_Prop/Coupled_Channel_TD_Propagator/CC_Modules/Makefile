#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
MAIN = $(MESA_HOME)
LIB = $(MAIN)/library
MOD_DIR = ../../Modules
PROP_MOD_DIR = ../Prop_Modules
PROP_MOD_LIB = $(PROP_MOD_DIR)/prop_mod_lib_f90.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c -O4  -I$(MOD_DIR)

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l
MOD_SRCS = \
	$(PROP_MOD_DIR)/auto_correlation_module.f90 \
	$(PROP_MOD_DIR)/moment_module.f90 \
	$(PROP_MOD_DIR)/plot_module.f90 \
	$(PROP_MOD_DIR)/regional_module.f90 \
	$(PROP_MOD_DIR)/normalize_module.f90 \
	$(PROP_MOD_DIR)/non_linear_potential_module.f90 \
	$(PROP_MOD_DIR)/initial_state_module.f90 \
	$(PROP_MOD_DIR)/check_convergence_module.f90 \
	$(PROP_MOD_DIR)/check_normalization_module.f90 \
	$(PROP_MOD_DIR)/test_overlap_module.f90 \
	$(PROP_MOD_DIR)/spatial_wavefunction_module.f90 \
	$(PROP_MOD_DIR)/matrix_vector_multiply_module.f90 \
	$(PROP_MOD_DIR)/arnoldi_module.f90 \
	$(PROP_MOD_DIR)/split_operator_module.f90 

%.mod : %.o
	@if [! -f $@ ]; then \
	rm $< \
	$(MAKE) $< \
	fi

%.o : %.f90
	$(FC) $(FFLAGS) $<

all: $(PROP_MOD_LIB)
$(PROP_MOD_LIB): $(MOD_SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(PROP_MOD_LIB) $(?:.f90=.o)
	$(RANLIB) $(PROP_MOD_LIB)

files: $(MOD_SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o *.mod *~
