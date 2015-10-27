#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
DVRLIB = ../dvrlib_f90.a
MODDIR =../../Modules
.SUFFIXES: .f90 $(SUFFIXES)
AR = ar
ARFLAGS = rvu

FFLAGS = -c -O4 -I$(MODDIR)
LD =$(FC)
LDFLAGS =
RCS = RCS
CO = co -l
REL = 

SRCS = \
	blcyl.f90 \
	blrad.f90 \
	blxyz.f90 \
	cnvtab.f90 \
	cnvtpt.f90 \
	genwtf.f90 \
	genq.f90 \
	genrwt.f90 \
	kemel.f90 \
	kinxyz.f90 \
	kinrad.f90 \
	kincyl.f90 \
	lancz.f90 \
	modab.f90 \
	plywpt.f90 \
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
	$(CO) $(RCS)/*.f90*

clean:
	rm *.o

