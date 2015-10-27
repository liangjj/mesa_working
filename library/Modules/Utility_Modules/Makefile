#    @(#)Makefile	5.1   11/6/94
#
#    Makefile for the utility portion of MesaLib
#    
include $(MODULES)/Makefile.inc
SRCS = \
	Matrix_Print.f90 
OBJECTS = \
	$(SRCS) $(SRCS:.f90=.o) 
#
$(UTIL_MOD_LIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $(UTIL_MOD_LIB) $(OBJECTS)
	$(RANLIB) $(UTIL_MOD_LIB)
	$(AR) $(ARFLAGS) $(MOD_LIB) $(OBJECTS)
	$(RANLIB) $(MOD_LIB)
sources:
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod
