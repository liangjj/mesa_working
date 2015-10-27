#@(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR = $(MESA_HOME)/Test_Array
include $(TOP_DIR)/Makefile.inc
SRCS = \
	Arrays.f90 

all: $(ARRAY_MODULE_LIB) 
$(ARRAY_MODULE_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(ARRAY_MODULE_LIB) $(?:.f90=.o)
	$(RANLIB) $(ARRAY_MODULE_LIB)



files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
#Make_Array.f90