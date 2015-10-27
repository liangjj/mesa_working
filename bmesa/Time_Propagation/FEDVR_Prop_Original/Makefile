#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
MAIN = /usr/home/bis/mesa
LIB = $(MAIN)/library
MODDIR = $(LIB)/Modules
MESABIN =  /usr/home/bis/mesa/bin
MODLIB = $(LIB)/modlib_f90.a
TIMPRP = $(MAIN)/bmesa/Time_Propagation
FEDVR_Prop = $(TIMPRP)/FEDVR_Prop
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
	add_absorb.f90 \
	auto_08_1d.f90 \
	auto_08_2d.f90 \
	auto_08_3d.f90 \
	chk_nrm_08_1d.f90 \
	chk_nrm_08_2d.f90 \
	chk_nrm_08_3d.f90 \
	cp_psi_08.f90 \
	c_vect_08.f90 \
	diag_reg_d.f90 \
	diag_reg_z.f90 \
	dvrprop_main.f90 \
	exp_off_diag.f90 \
	exp_off_diag_m_v.f90 \
	exp_off_diag_v_m_2_d.f90 \
	exp_off_diag_v_m_3_d.f90 \
	fil_reg_dvr.f90 \
	fil_reg_fd.f90 \
	fil_08.f90 \
	gauss_paket_08.f90 \
	ke_reg.f90 \
	modify_diag.f90 \
	moment_08.f90 \
	mom_1d.f90 \
	mom_2d.f90 \
	mom_3d.f90 \
	plot_08.f90 \
	plot_1d.f90 \
	plot_2d.f90 \
	plot_3d.f90 \
        plot_propagator.f90 \
        plot_wavefunction.f90 \
	pert.f90 \
	propagator_2.f90 \
	rad_paket_08.f90 \
	real_space_propagator_2_order.f90 \
	so_prop.f90 \
	sppose_08.f90 \
	umat_reg_d.f90 \
	umat_reg_z.f90 \
	v_m_v_gen.f90 \
	v_v_m_gen.f90 \
	v_m_v_2.f90 \
	v_m_v_3.f90 \
	v_m_v_4.f90 \
	v_m_v_5.f90 \
	v_m_v_6.f90 \
	v_m_v_7.f90 \
	v_m_v_8.f90 \
	v_m_v_9.f90 \
	v_m_v_10.f90 \
	v_v_m_2.f90 \
	v_v_m_3.f90 \
	v_v_m_4.f90 \
	v_v_m_5.f90 \
	v_v_m_6.f90 \
	v_v_m_7.f90 \
	v_v_m_8.f90 \
	v_v_m_9.f90 \
	v_v_m_10.f90 \
	z_fil_08.f90 \
	z_proj_08.f90

.f90.o:
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


sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o

