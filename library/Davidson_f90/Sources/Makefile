#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_LIB)/Davidson_f90
include $(TOP_DIR)/Makefile.inc
SRCS = \
	Read_Matrix_From_Input.f90



.mod.o:

.f90.o:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(DAVIDSON_MODULES) $<
.f90.mod:
	$(F90) $(F90FLAGS) -I$(MODULES) -I$(DAVIDSON_MODULES) $<

all: $(DAVIDSON_SOURCE_LIB)
$(DAVIDSON_SOURCE_LIB): $(SRCS:.f90=.o)   
	$(AR) $(ARFLAGS) $(DAVIDSON_SOURCE_LIB) $(?:.f90=.o)
	$(RANLIB) $(DAVIDSON_SOURCE_LIB)

files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *~
