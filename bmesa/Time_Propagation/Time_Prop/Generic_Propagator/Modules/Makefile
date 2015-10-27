#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#

.SUFFIXES: .f90 .mod $(SUFFIXES)
GENERIC_PROPAGATOR_DIR = \
$(MESA_HOME)/bmesa/Time_Propagation/Time_Prop/Generic_Propagator
MOD_DIR = $(GENERIC_PROPAGATOR_DIR)/Modules
MOD_LIB = $(MOD_DIR)/mod_lib_f90.a

RCS = RCS
CO = co -l

SRCS = \
	input_output.f90 \
	potential.f90 \
	grid_global.f90 \
	dvr_shared.f90 \
	dvr_global.f90 \
	prop_prnt.f90 \
	prop_global.f90 \
	arnoldi_global.f90\
	arnoldi_global_it.f90 \
	arnoldi_global_rt.f90 \
	dvd_prnt.f90 \
	dvr_prnt.f90 \
	dvd_global.f90 \
	fd_prnt.f90 \
	fd_global.f90 \
	lanczos_prnt.f90 \
	lanczos_global.f90 \
	dvrprop_global_rt.f90 \
	dvrprop_global_it.f90 \
	dvrprop_global.f90 \
	r_matrix_global.f90 \
	mesa_global.f90

.f90.o:
	$(FC) $(FFLAGS) -I$(MOD_DIR) $<

all: $(MOD_LIB) 
$(MOD_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(MOD_LIB) $(?:.f90=.o)
sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.a

