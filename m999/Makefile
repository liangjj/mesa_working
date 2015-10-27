#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
MAIN = $(MESA_HOME)
LIB = $(MAIN)/library
LIBMODDIR = $(LIB)/Modules
MESABIN =  $(MESA_BIN)
MODLIB = $(LIB)/modlib_f90.a
MESALIB = $(LIB)/mesalib.a

FFLAGS = -c -O5 -I$(LIBMODDIR) 
LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	test_type.f90 \
	test_mul.f90

.f90.o:
	$(FC) $(FFLAGS) $<

all: $(MESABIN)/m999
$(MESABIN)/m999: $(SRCS:.f90=.o) \
                  $(MODLIB) $(MESALIB) $(MDLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f90=.o) \
        $(MODLIB) $(MESALIB) $(MDLIB) -o $(MESABIN)/m999     



