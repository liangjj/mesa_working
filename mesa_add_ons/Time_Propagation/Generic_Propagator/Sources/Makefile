#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/Generic_Propagator
include $(TOP_DIR)/Makefile.inc
SRCS = \
	pt_reg.f90 \
	add_absorb.f90 \
	ke_reg_dvr.f90 \
	ke_reg_fd.f90 \
	input_prop.f90 \
	space_prop.f90 \
	modify_diag.f90 \
	check_gs_energy.f90 \
	v_couple.f90 \
	v_tim.f90 \
	mk_phi.f90 \
	nr_paket.f90 \
	ov1_quad.f90 \
	ov2_quad.f90 \
	pert.f90 \
	sc_dip.f90 

RCS = RCS
CO = co -l

.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(PROP_MODULES) $<

all: $(PROP_SOURCE_LIB)
$(PROP_SOURCE_LIB): $(SRCS:.f90=.o)   
	$(AR) $(ARFLAGS) $(PROP_SOURCE_LIB) $(?:.f90=.o)
	$(RANLIB) $(PROP_SOURCE_LIB)

files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o
	rm -f *.a
