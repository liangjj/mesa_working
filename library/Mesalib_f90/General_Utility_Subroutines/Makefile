#    @(#)Makefile	5.1   11/6/94
#
#    Makefile for the utility portion of MesaLib
#    
TOP_DIR = $(MESA_LIB)/Mesalib_f90
include $(TOP_DIR)/Makefile.inc
SRCS = \
	Matrix_Print.f90 

all: $(UTIL_LIB) 
$(UTIL_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(UTIL_LIB) $(?:.f90=.o)
	$(AR) $(ARFLAGS) $(UTILITY_LIB) $(?:.f90=.o)

sources:
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o
