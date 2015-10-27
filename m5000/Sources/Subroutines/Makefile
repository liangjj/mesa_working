#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/m5000
include $(TOP_DIR)/Makefile.inc
include $(MESA_HOME)/Makefile.inc
SRCS = \
	pt_reg.f90 \
	add_absorb.f90 \
	modify_diag.f90 \
	check_gs_energy.f90 \
	v_couple.f90 \
	v_tim.f90 \
	ke_reg_dvr.f90 \
	ke_reg_fd.f90 \
	mk_phi.f90 \
	nr_paket.f90 \
	ov1_quad.f90 \
	ov2_quad.f90 \
	pert.f90 \
	sc_dip.f90 \
	Laser_Field.f90
RCS = RCS
CO = co -l

all: $(SIL_SUBROUTINES_LIB)
$(SIL_SUBROUTINES_LIB): $(SRCS:.f90=.o)   
	$(AR) $(ARFLAGS) $(SIL_SUBROUTINES_LIB) $(?:.f90=.o)
	$(RANLIB) $(SIL_SUBROUTINES_LIB)
clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(OBJECT_DIR)/*.o
	rm -f $(MODULE_DIR)/*.mod
	rm -f $(LIB_DIR)/*.a
movein:
	mv *.mod $(MODULE_DIR)/.
	mv *.o $(OBJECT_DIR)/.

moveout:
	mv $(MODULE_DIR)/*.mod .
	mv $(OBJECT_DIR)/*.o .
