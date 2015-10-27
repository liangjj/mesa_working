#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
MAIN = $(MESA_HOME)
LIB = $(MAIN)/library
MODDIR = $(LIB)/Modules
TIMPRP = $(MAIN)/bmesa/Time_Propagation
IT_MODULES = $(TIMPRP)/it_modules_f90.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c -O5 

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l
SRCS = \
	auto_correlation.f90 \
	gaussian_wave_packet.f90 \
	h_on_vector.f90 \
	plot_wavefunction.f90 \
	normalize.f90 \
	eigenvalues.f90 \
	moment.f90 \
	psi_h_psi.f90 \
	initial_state.f90 \
	non_linear_potential.f90 \
	v_m_v.f90 \
	v_v_m.f90 \
	exp_off_diagonal.f90 \
	real_space_propagator.f90 \
	finite_element_propagator.f90 \
	finite_element_matrix_multiply.f90 \
	diagonal_mul.f90 \
	diag_reg.f90 \
	modify_diag.f90 \
        plot_propagator.f90 \
	propagator.f90 \
	pt_reg.f90 \
	setup.f90 \
	umat_reg.f90 \
	z_proj.f90 \
	exp_diag_prop.f90 \
	ke_reg_dvr.f90 \
	ke_reg_fd.f90 \
	dvrprop_main_it.f90 

.f90.o:
	$(FC) $(FFLAGS) -I$(MODDIR) $<

all: $(IT_MODULES)
$(IT_MODULES): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(IT_MODULES) $(?:.f90=.o)     

sources: $(SRCS) $(LIBSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o *.mod *~
