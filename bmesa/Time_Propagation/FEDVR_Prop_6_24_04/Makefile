#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
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

MODSRCS = \
	auto_correlation.f90 \
	moment.f90 \
	normalize.f90 \
	plot.f90 \
	regional_diag.f90 \
	regional_umat.f90 \
	v_m_v.f90 \
	v_v_m.f90 \
	exp_off_diagonal.f90 
SRCS = \
	add_absorb.f90 \
	cp_psi.f90 \
	c_vect.f90 \
	diag_reg_d.f90 \
	diag_reg_z.f90 \
	dvrprop_main.f90 \
	fill_d.f90 \
	gauss_paket.f90 \
	ke_reg_dvr.f90 \
	ke_reg_fd.f90 \
	modify_diag.f90 \
        plot_propagator.f90 \
	pert.f90 \
	propagator_2.f90 \
	propagator_4.f90 \
	rad_paket.f90 \
	real_space_propagator_2_order.f90 \
	so_prop.f90 \
	sppose.f90 \
	umat_reg_d.f90 \
	umat_reg_z.f90 \
	z_fill.f90 \
	z_proj.f90

%.mod : %.o
	@if [! -f $@ ]; then \
	rm $< \
	$(MAKE) $< \
        fi

%.o : %.f90
	$(FC) $(FFLAGS) -I$(MODDIR) $<

all: $(MESABIN)/m1104
$(MESABIN)/m1104: $(MODSRCS:.f90=.o) $(SRCS:.f90=.o)  \
                  $(TIMLIB) $(TPRPLIB) $(DVRLIB) $(FDLIB) \
                  $(MODLIB) $(HAMLIB) $(POTLIB) $(ITDVRLIB) \
                  $(TDVRLIB) $(MESALIB) $(MDLIB)
	$(LD) $(LDFLAGS) $(MODSRCS:.f90=.o) $(SRCS:.f90=.o)  \
	$(TIMLIB) $(TPRPLIB) $(DVRLIB) $(FDLIB) $(MODLIB) \
	$(HAMLIB) $(POTLIB) $(ITDVRLIB) $(TDVRLIB) \
	$(MESALIB) $(MDLIB) -o $(MESABIN)/m1104     

sources: $(SRCS) $(LIBSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o *.mod *~
