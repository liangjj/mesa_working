#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
MAIN = $(MESA_HOME)
LIB = $(MAIN)/library
MODDIR = $(LIB)/Modules
MESABIN =  $(MESA_BIN)
MODLIB = $(LIB)/modlib_f90.a
TIMPRP = $(MAIN)/bmesa/Time_Propagation
DVRLIB = $(LIB)/dvrlib_f90/dvrlib_f90.a
FDLIB = $(LIB)/fdlib_f90/fdlib_f90.a
POTLIB = $(LIB)/potlib/potlib_f90.a
HAMLIB = $(LIB)/hamlib/hamlib.a
MESALIB = $(LIB)/mesalib.a
ITDVRLIB  = $(LIB)/itlib/itdvrlib.a
TDVRLIB = $(LIB)/tdvrlib/tdvrlib.a
TIMLIB = $(TIMPRP)/timlib_f90.a
MDLIB = /usr/lib/libcxml.a    

FC = fort
FFLAGS = -c -O5 -I$(MODDIR)

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

ARNSRCS = \
        check_convergence.f90 \
        check_normalization.f90 \
        diagonalize.f90 \
        exponential_propagator.f90 \
        small_matrices.f90 \
        test_overlap.f90 \
        trial_vectors.f90


MODSRCS = \
	auto_correlation.f90 \
	gaussian_wave_packet.f90 \
	non_linear_potential.f90 \
	plot_wavefunction.f90 \
	normalize.f90 \
	eigenvalues.f90 \
	moment.f90 \
	psi_h_psi.f90 \
	initial_state.f90 \
	v_out_eq_dvr_mat_v_in.f90 \
	v_out_eq_v_in_dvr_mat.f90 \
	exp_off_diagonal.f90 \
	real_space_propagator.f90 \
	pack_and_read_h.f90 \
	packed_dvr_matrix_vector_multiply.f90 \
	dvr_matrix_vector_multiply.f90 \
	finite_element_matrix_multiply.f90 \
	fd_m_v.f90 \
	h_m_v.f90 \
	h_on_vector.f90 \
        initialize.f90 \
        arnoldi_driver.f90 \
	finite_element_propagator.f90 



SRCS = \
	arnoldi_dat.f90 \
	arnoldi_propagator.f90 \
	diagonal_mul.f90 \
	diag_reg.f90 \
	exp_diag_prop.f90 \
	input_prop.f90 \
	ke_reg_dvr.f90 \
	ke_reg_fd.f90 \
	modify_diag.f90 \
        plot_propagator.f90 \
	propagator.f90 \
	pt_reg.f90 \
	regional_matrices.f90 \
	setup.f90 \
	space_prop.f90 \
	split_operator.f90 \
	time_prop.f90 \
	umat_reg.f90 \
	z_proj.f90 \
	prop_main_it.f90

%.mod : %.o
	@if [! -f $@ ]; then \
	rm $< \
	$(MAKE) $< \
        fi

%.o : %.f90
	$(FC) $(FFLAGS) $<

all: $(MESABIN)/m1200
$(MESABIN)/m1200: $(ARNSRCS:.f90=.o)           \
                  $(MODSRCS:.f90=.o) \
		  $(SRCS:.f90=.o)  \
                  $(TIMLIB) $(DVRLIB) $(FDLIB) \
                  $(MODLIB) $(HAMLIB) $(POTLIB) $(ITDVRLIB) \
                  $(TDVRLIB) $(MESALIB) $(MDLIB)
	$(LD) $(LDFLAGS) $(ARNSRCS:.f90=.o)          \
                         $(MODSRCS:.f90=.o) \
			 $(SRCS:.f90=.o)  \
	$(TIMLIB) $(DVRLIB) $(FDLIB) $(MODLIB) \
	$(HAMLIB) $(POTLIB) $(ITDVRLIB) $(TDVRLIB) \
	$(MESALIB) $(MDLIB) -o $(MESABIN)/m1200     

sources: $(SRCS) $(LIBSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o *.mod *~
