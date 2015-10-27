#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
SOURCE = m5001
include $(MESA_HOME)/$(SOURCE)/Makefile.inc
SRCLOC = $(MESA_HOME)/$(SOURCE)/Subroutines
SRCS = \
	$(SRCLOC)/pt_reg.f90 \
	$(SRCLOC)/add_absorb.f90 \
	$(SRCLOC)/modify_diag.f90 \
	$(SRCLOC)/check_gs_energy.f90 \
	$(SRCLOC)/v_couple.f90 \
	$(SRCLOC)/v_tim.f90 \
	$(SRCLOC)/ke_reg_dvr.f90 \
	$(SRCLOC)/ke_reg_fd.f90 \
	$(SRCLOC)/mk_phi.f90 \
	$(SRCLOC)/nr_paket.f90 \
	$(SRCLOC)/ov1_quad.f90 \
	$(SRCLOC)/ov2_quad.f90 \
	$(SRCLOC)/pert.f90 \
	$(SRCLOC)/sc_dip.f90 \
	$(SRCLOC)/Laser_Field.f90

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
