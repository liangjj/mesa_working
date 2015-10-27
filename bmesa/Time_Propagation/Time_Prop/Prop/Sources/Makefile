#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
MAIN = $(MESA_HOME)
LIB = $(MAIN)/library
LIBMODDIR = $(LIB)/Modules
MODDIR = ../Modules
SOURCELIB = ../Sources/sourcelib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

SRCS = \
	../Sources/pt_reg.f90 \
	../Sources/add_absorb.f90 \
	../Sources/ke_reg_dvr.f90 \
	../Sources/ke_reg_fd.f90 \
	../Sources/input_prop.f90 \
	../Sources/space_prop.f90 \
	../Sources/modify_diag.f90 \
	../Sources/setup.f90 \
	../Sources/ground_state_energy.f90 \
	../Sources/check_gs_energy.f90 \
	../Sources/v_couple.f90 \
	../Sources/v_tim.f90 \
	../Sources/mk_phi.f90 \
	../Sources/nr_paket.f90 \
	../Sources/ov1_quad.f90 \
	../Sources/ov2_quad.f90 \
	../Sources/pert.f90 \
	../Sources/sc_dip.f90 

FFLAGS = -c -O5 -I$(LIBMODDIR) -I$(MODDIR)

RCS = RCS
CO = co -l

.f90.o:
	$(FC) $(FFLAGS) $<

all: $(SOURCELIB)
$(SOURCELIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(SOURCELIB) $(?:.f90=.o)
	$(RANLIB) $(SOURCELIB)

files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm *.o *.mod *~
