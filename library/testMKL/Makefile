#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
MAIN=/redhat/home/mesa/mesa/library/testMKL
SRCS =	Test.f90 
LIB_PATH =/opt/intel/mkl/10.0.1.014/lib/em64t
F90FLAGS = -c -O4
.f90.o:
	$(F90) $(F90FLAGS) $<

all: $(MAIN)/Test
$(MAIN)/Test: $(SRCS:.f90=.o) $(LIB_PATH) $(MD_LIB)
	$(LD90) $(LDFLAGS) -L$(LIB_PATH) $(SRCS:.f90=.o) -l$(MD_LIB) -o $(MAIN)/Test

files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
