#   @(#)Makefile	1.4   8/2/91
#
# Makefile for TPRPLIB
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
FFLAGS = -c -O4 -I$(MODDIR)
MAIN = $(MESA_HOME)
LIB = $(MAIN)/library
MODDIR = $(LIB)/Modules
TIMPRP = $(MAIN)/bmesa/Time_Propagation
TIMLIB = $(TIMPRP)/timlib_f90.a

GET = sccs get

SRCS = \
	cp_psi0.f90 \
	c_vect0.f90 \
	fil_1.f90 \
	fil_2.f90 \
	fil_3.f90 \
	gauss_paket.f90 \
	moment.f90 \
	mk_phi.f90 \
	ov1_quad.f90 \
	ov2_quad.f90 \
	nr_paket.f90 \
	pack_h.f90 \
	pert.f90 \
	psi0_fil_2.f90 \
	psi0_fil_3.f90 \
	rad_paket.f90 \
	sc_dip.f90 \
	soln.f90 \
	sppose.f90 \
	trial.f90 \
	v_couple.f90 \
	v_dw_plus.f90 \
	v_dw_minus.f90 \
	v_nl_rt_1d.f90 \
	v_nl_it_1d.f90 \
	v_nl_rt_2d.f90 \
	v_nl_it_2d.f90 \
	v_nl_rt_3d.f90 \
	v_nl_it_3d.f90 \
	vnlse.f90 \
	v_tim.f90 \
	zfil_2.f90 \
	zfil_3.f90

.f90.o:
	$(FC) $(FFLAGS) $<

all: $(TIMLIB)
$(TIMLIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(TIMLIB) $(?:.f90=.o)
