#    @(#)Makefile	5.2   4/17/95
#
#    Makefile for the mathematics library of MesaLib
#    
TOP_DIR = $(MESA_LIB)/Mesalib_f90
include $(TOP_DIR)/Makefile.inc
#
SRCS = \
	Associated_Legendre_Functions.f90 

all: $(MATH_LIB)
$(MATH_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(MATH_LIB) $(?:.f90=.o)
	$(AR) $(ARFLAGS) $(UTILITY_LIB) $(?:.f90=.o)

sources: $(SRCS)
$(SRCS):
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o

