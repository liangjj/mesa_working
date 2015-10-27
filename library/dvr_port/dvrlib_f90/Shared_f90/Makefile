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
	blcyl.f90 \
	blrad.f90 \
	blxyz.f90 \
	cnvtab.f90 \
	cnvtpt.f90 \
	diag_reg.f90 \
	genq.f90 \
	genrwt.f90 \
	genwtf.f90 \
	kemel.f90 \
	ke_reg_dvr.f90 \
	ke_reg_fd.f90 \
	kincyl.f90 \
	kinrad.f90 \
	kinxyz.f90 \
	lancz.f90 \
	modab.f90 \
	modify_diag.f90 \
	pt_reg.f90 \
	recur.f90 \
	scaprd.f90 \
	todelta.f90 \
	xsq2x.f90
.f90.o:
	$(FC) $(FFLAGS) $<

all: $(DVRLIB)
$(DVRLIB): $(SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(DVRLIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o

