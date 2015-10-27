#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator input
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
MAIN = $(MESA_HOME)
LIB = $(MAIN)/library
MODDIR = $(LIB)/Modules
MESABIN =  $(MESA_BIN)
MODLIB = $(LIB)/modlib_f90.a
TIMPRP = $(MAIN)/bmesa/Time_Propagation
FEDVR_Prop_RT = $(TIMPRP)/FEDVR_Prop_RT
DVRLIB = $(LIB)/dvrlib_f90/dvrlib_f90.a
FDLIB = $(LIB)/fdlib_f90/fdlib_f90.a
POTLIB = $(LIB)/potlib/potlib_f90.a
HAMLIB = $(LIB)/hamlib/hamlib.a
MESALIB = $(LIB)/mesalib.a
ITDVRLIB  = $(LIB)/itlib/itdvrlib.a
TDVRLIB = $(LIB)/tdvrlib/tdvrlib.a
TIMLIB = $(TIMPRP)/timlib_f90.a
TPRPLIB = $(TIMPRP)/tprplib_f90.a
MDLIB = /usr/lib/libcxml.a    

FC = fort
FFLAGS = -c -O5 

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	auto_correlation.f90 \
	plot_wavefunction.f90 \
	gaussian_wave_packet.f90 \
	normalize.f90 \
	moment.f90 \
	initial_state.f90 \
	non_linear_potential.f90 \
	regional_diag.f90 \
	regional_umat.f90 \
	v_m_v_rt.f90 \
	v_v_m_rt.f90 \
	exp_off_diagonal.f90 \
	real_space_propagator.f90 \
	finite_element_propagator.f90 \
	add_absorb.f90 \
	cs_diag.f90 \
	diagonal_mul.f90 \
	propagator.f90 \
	dvrprop_main_rt.f90 \
	evolve.f90 \
	ke_reg_dvr.f90 \
	ke_reg_fd.f90 \
	modify_diag.f90 \
        plot_propagator.f90 \
	pt_reg.f90 \
	setup.f90 \
	z_proj.f90 

%.mod : %.o
	@if [! -f $@ ]; then \
	rm $< \
	$(MAKE) $< \
        fi

%.o : %.f90
	$(FC) $(FFLAGS) -I$(MODDIR) $<

all: $(MESABIN)/m1104
$(MESABIN)/m1104: $(SRCS:.f90=.o)  \
                  $(TIMLIB) $(TPRPLIB) $(DVRLIB) $(FDLIB) \
                  $(MODLIB) $(HAMLIB) $(POTLIB) $(ITDVRLIB) \
                  $(TDVRLIB) $(MESALIB) $(MDLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f90=.o)  \
	$(TIMLIB) $(TPRPLIB) $(DVRLIB) $(FDLIB) $(MODLIB) \
	$(HAMLIB) $(POTLIB) $(ITDVRLIB) $(TDVRLIB) \
	$(MESALIB) $(MDLIB) -o $(MESABIN)/m1104     

sources: $(SRCS) $(LIBSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o *.mod *~
