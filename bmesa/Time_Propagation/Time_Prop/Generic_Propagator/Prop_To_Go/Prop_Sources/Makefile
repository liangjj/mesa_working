#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
GENERIC_PROPAGATOR_DIR = \
/home/mesa/mesa/bmesa/Time_Propagation/Time_Prop/Generic_Propagator
MOD_DIR = $(GENERIC_PROPAGATOR_DIR)/Modules
MOD_LIB = $(MOD_DIR)/mod_lib_f90.a
PROP_MOD_DIR = $(GENERIC_PROPAGATOR_DIR)/Prop_To_Go/Prop_Modules
PROP_MOD_LIB = $(PROP_MOD_DIR)/prop_mod_lib_f90.a
PROP_SOURCES_DIR = $(GENERIC_PROPAGATOR_DIR)/Prop_To_Go/Prop_Sources
PROP_SOURCE_LIB = $(PROP_SOURCES_DIR)/prop_source_lib_f90.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

SRCS = \
	$(PROP_SOURCES_DIR)/pt_reg.f90 \
	$(PROP_SOURCES_DIR)/add_absorb.f90 \
	$(PROP_SOURCES_DIR)/ke_reg_dvr.f90 \
	$(PROP_SOURCES_DIR)/ke_reg_fd.f90 \
	$(PROP_SOURCES_DIR)/input_prop.f90 \
	$(PROP_SOURCES_DIR)/space_prop.f90 \
	$(PROP_SOURCES_DIR)/modify_diag.f90 \
	$(PROP_SOURCES_DIR)/check_gs_energy.f90 \
	$(PROP_SOURCES_DIR)/v_couple.f90 \
	$(PROP_SOURCES_DIR)/v_tim.f90 \
	$(PROP_SOURCES_DIR)/mk_phi.f90 \
	$(PROP_SOURCES_DIR)/nr_paket.f90 \
	$(PROP_SOURCES_DIR)/ov1_quad.f90 \
	$(PROP_SOURCES_DIR)/ov2_quad.f90 \
	$(PROP_SOURCES_DIR)/pert.f90 \
	$(PROP_SOURCES_DIR)/sc_dip.f90 

RCS = RCS
CO = co -l

.f90.o:
	$(FC) $(FFLAGS) -I$(MOD_DIR) -I$(PROP_MOD_DIR) $<

all: $(PROP_SOURCE_LIB)
$(PROP_SOURCE_LIB): $(SRCS:.f90=.o)   
	$(AR) $(ARFLAGS) $(PROP_SOURCE_LIB) $(?:.f90=.o)
	$(RANLIB) $(PROP_SOURCE_LIB)

files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o
	rm -f *.a
