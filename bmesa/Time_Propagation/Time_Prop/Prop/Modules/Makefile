#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
MAIN = $(MESA_HOME)
LIB = $(MAIN)/library
LIBMODDIR = $(LIB)/Modules
PROPMODLIB = ../Modules/propmodlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c -O5 -I$(LIBMODDIR)

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

MODSRCS = \
	../Modules/moment_module.f90 \
	../Modules/regional_module.f90 \
	../Modules/auto_correlation_module.f90 \
	../Modules/normalize_module.f90 \
	../Modules/plot_module.f90 \
	../Modules/non_linear_potential_module.f90 \
	../Modules/initial_state_module.f90 \
	../Modules/check_convergence_module.f90 \
	../Modules/check_normalization_module.f90 \
	../Modules/test_overlap_module.f90 \
	../Modules/spatial_wavefunction_module.f90 \
	../Modules/matrix_vector_multiply_module.f90 \
	../Modules/arnoldi_module.f90 \
	../Modules/split_operator_module.f90 

%.mod : %.o
	@if [! -f $@ ]; then \
	rm $< \
	$(MAKE) $< \
	fi

%.o : %.f90
	$(FC) $(FFLAGS) $<

all: $(PROPMODLIB)
$(PROPMODLIB): $(MODSRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(PROPMODLIB) $(?:.f90=.o)
	$(RANLIB) $(PROPMODLIB)

files: $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o *.mod *~
