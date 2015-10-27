#   @(#)Makefile	1.4   8/2/91
#
#
include ../Makefile.inc

SRCS = \
	Utility_Module.f90 \
	Read_DVR_Module.f90 \
	Renormalization_Module.f90 \
	DVR_Polynomials_Module.f90 \
	DVR_Kinetic_Energy_Module.f90 \
	DVR_H_0_Module.f90 \
	DVR_Module.f90 \
	Matrix_Scale_and_Assemble.f90 \
	Matrix_Diagonalization.f90 \
	Poisson_Module.f90 \
	P_LM.f90 \
	Two_Electron_FEDVR_Module.f90

all: $(FEDVR_MODULE_LIB)
$(FEDVR_MODULE_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(FEDVR_MODULE_LIB) $(?:.f90=.o)
	$(RANLIB) $(FEDVR_MODULE_LIB)
	mv *.mod $(MODULE_DIRECTORY)/.
	mv *.o $(OBJECT_DIRECTORY)/.

clean:
	rm -f *.o 
	rm -f *.mod 
	rm -f $(OBJECT_DIRECTORY)/*.o 
	rm -f $(MODULE_DIRECTORY)/*.mod
	rm -f $(LIB_DIRECTORY)/*.a
movein:
	mv *.mod $(MODULE_DIRECTORY)/.
	mv *.o $(OBJECT_DIRECTORY)/.

moveout:
	mv $(MODULE_DIRECTORY)/*.mod .
	mv $(OBJECT_DIRECTORY)/*.o .