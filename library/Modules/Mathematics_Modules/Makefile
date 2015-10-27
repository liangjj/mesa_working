#    @(#)Makefile	5.2   4/17/95
#
#    Makefile for the mathematics library of MesaLib
#    
include $(MODULES)/Makefile.inc
.SUFFIXES: .f90 .mod $(SUFFIXES)
#
SRCS = \
	Special_Functions.f90 \
	Lentz_Thompson.f90 \
	Prolate_Functions.f90 \
	Associated_Legendre_Functions.f90 \
	Regular_Associated_Legendre_Functions.f90 
OBJECTS = \
	$(SRCS:.f90=.o) 
#
$(MATH_MOD_LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $(MATH_MOD_LIB) $(OBJECTS)
	$(RANLIB) $(MATH_MOD_LIB)
	$(AR) $(ARFLAGS) $(MOD_LIB) $(OBJECTS)
	$(RANLIB) $(MOD_LIB)
sources: $(SRCS)

clean:
	rm -f *.o *.mod
