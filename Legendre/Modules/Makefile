#@(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR = $(MESA_HOME)/Legendre
include $(TOP_DIR)/Makefile.inc
SRCS = \
	itoc.f \
	input_output.f90 \
	Matrix_Print.f90 \
	accuracy.f90 \
	Data_Module.f90 \
	Special_Functions.f90 \
	Lentz_Thompson.f90 \
	Associated_Legendre_Functions.f90 

all: $(LEGENDRE_MODULES_LIB) 
$(LEGENDRE_MODULES_LIB): $(SRCS:.f90=.o) $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(LEGENDRE_MODULES_LIB) $(?:.f90=.o)
	$(RANLIB) $(LEGENDRE_MODULES_LIB)



files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
#	Associated_Legendre_Functions.f90 \