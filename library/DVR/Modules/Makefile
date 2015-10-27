#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
.SUFFIXES: .f90 $(SUFFIXES)
include $(TOP_DIR)/Makefile.inc
RCS = RCS
CO = co -l
REL = 

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
	fourier.f90 \
	fourier_basis.f90 \
	hamphy.f90 \
	kemat.f90 \
	ke_spherical.f90 \
	ke_hermite.f90 \
	ke_laguerre.f90 \
	ke_legendre.f90 \
	ke_cyl.f90 \
	ke_theta.f90 \
	lgngr.f90 \
	lobatto.f90 \
	makd.f90 \
	makd_1.f90 \
	makod.f90 \
	makod_1.f90 \
	memreg.f90 \
	ovmat.f90 \
	p_hermite.f90 \
	phy.f90 \
	p_laguerre.f90 \
	pmat.f90 \
	poisson.f90 \
	potntl.f90 \
	rhofun.f90 \
	rmamp.f90 \
	twomat.f90 \
	vrmat.f90


all: $(SPHEROIDAL_DVR_LIB)
$(SPHEROIDAL_DVR_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(SPHEROIDAL_DVR_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90*

clean:
	rm -f *.o
