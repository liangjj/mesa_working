# examples/f90/Makefile.  Generated from Makefile.in by configure.
SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .o .f90 .txt 
PROG = pi3f90
F90 = mpif90
SRCS =  pi3f90.f90
OBJS =  $(SRCS:.f90=.o)
LIBS = -lpthread -lrt -lmpich  
.f90.o:
	$(F90) $(F90FLAGS) $<
.f.o:
	$(F77) $(FFLAGS) $<
all: $(PROG)
$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) $(OBJS) $(LIBS) $(UTILITY_LIB) $(MD_LIB) -o $@
clean:
	-rm -f $(OBJS)


