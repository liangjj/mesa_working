#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
BIN = ../Matrix_Exponential
FFLAGS = -c -O5 
LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	input_output.f90 \
	prntcmn.f90 \
	prntfmn.f90 \
	ebct.f90 \
	cebct.f90 \
	Matrix_Exponential_Module.f90 \
	Matrix.f90

.f90.o:
	$(FC) $(FFLAGS) $<

all: $(BIN)/mat_exp
$(BIN)/mat_exp: $(SRCS:.f90=.o) $(MDLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f90=.o) $(MDLIB) \
	-o $(BIN)/mat_exp     



