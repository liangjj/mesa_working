#@(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
#  Replace the next two equalities with your fortran compiler
#  and compiler flags.

LEGENDRE_DRIVER = ./legendre_driver_x 

NOLINK = -c
.f90.o:
	$(F90) $(NOLINK) $(F90FLAGS) $(MODULE_LIST) $(F90FLAGS) $<
LIB = $(MOD_LIB) $(POTENTIAL_LIB) $(UTILITY_LIB)

LEG_DRIVER = \
	Legendre_Function_Driver.f90 

all: $(LEGENDRE_DRIVER) 
$(LEGENDRE_DRIVER): $(SRCS:.f90=.o) $(LEG_DRIVER:.f90=.o) 
	$(F90) $(F90FLAGS) $(MODULE_LIST) $(SRCS:.f90=.o) $(LEG_DRIVER:.f90=.o) $(LIB) -o $(LEGENDRE_DRIVER)

files: $(SRCS) $(LEG_DRIVER)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
