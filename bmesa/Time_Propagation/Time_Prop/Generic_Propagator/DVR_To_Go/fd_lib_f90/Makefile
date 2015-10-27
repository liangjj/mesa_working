#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
FD_LIB = ../fd_lib_f90.a
PROP_MOD_DIR = ../../Prop_To_Go/Modules
MOD_DIR = ../../Modules
AR = ar
ARFLAGS = rvu
RANLIB = ranlib

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	band3.f90 \
        band5.f90 \
	band7.f90 \
	cpy_3.f90 \
	fd_input.f90 \
	fd_basis.f90 \
	points.f90 

.f90.o:
	$(FC) $(FFLAGS) -I$(MOD_DIR) $<

all: $(FD_LIB)
$(FD_LIB): $(SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(FD_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.a

