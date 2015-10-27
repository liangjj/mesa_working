#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#

DVR_PORT = $(MESA_LIB)/dvr_port
DVRLIB = $(DVR_PORT)/dvrlib_f90.a
MODDIR = $(DVR_PORT)/Modules

.SUFFIXES: .f90 .mod $(SUFFIXES)

AR = ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c -O4 -I$(MODDIR)

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	blmat.f90 \
	conbsis.f90 \
	conham.f90 \
	conke.f90 \
	conv.f90 \
	cordfn.f90 \
	cpoly.f90 \
	diarep.f90 \
	drvply.f90 \
	filfun.f90 \
	filv.f90 \
	fourier_basis.f90 \
	fourier.f90 \
	hamphy.f90 \
	kemat.f90 \
	ke_hermite.f90 \
	ke_laguerre.f90 \
	ke_legendre.f90 \
	lgngr.f90 \
	lobatto.f90 \
	makd_1.f90 \
	makd.f90 \
	makod_1.f90 \
	makod.f90 \
	memreg.f90 \
	ovmat.f90 \
	p_hermite.f90 \
	phy.f90 \
	p_laguerre.f90 \
	pmat.f90 \
	poisson.f90 \
	rhofun.f90 \
	rmamp.f90 \
	twomat.f90 \
	vrmat.f90
.f90.o:
	$(FC) $(FFLAGS) $<

all: $(DVRLIB)
$(DVRLIB): $(SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(DVRLIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o

