#   @(#)Makefile	1.4   8/2/91
#
# Makefile for TIMLIB
#
FDLIB =fdlib_f90.a
MODDIR = ../Modules
.SUFFIXES: .f90 .mod $(SUFFIXES)
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
FFLAGS = -c -O4 -I$(MODDIR)

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
	$(FC) $(FFLAGS) -I$(MODDIR) $<

all: $(FDLIB)
$(FDLIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(FDLIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90
