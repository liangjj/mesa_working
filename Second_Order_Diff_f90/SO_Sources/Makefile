#   @(#)Makefile	1.4   8/2/91
#
TOP_DIR = $(MESA_HOME)/Second_Order_Diff_f90
include $(TOP_DIR)/Makefile.inc

SRCS = \
	add_absorb.f90 \
	ke_reg_dvr.f90 \
	ke_reg_fd.f90 \
	Make_Hamiltonian.f90 \
	modify_diag.f90 \
	pt_reg.f90 \
	Phase_Shift.f90 \
	Space_DVR.f90 \
	v_mat.f90

RCS = RCS
CO = co -l

all: $(SO_SOURCE_LIB)
$(SO_SOURCE_LIB): $(SRCS:.f90=.o)   
	$(AR) $(ARFLAGS) $(SO_SOURCE_LIB) $(?:.f90=.o)
	$(RANLIB) $(SO_SOURCE_LIB)

files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm *.o
