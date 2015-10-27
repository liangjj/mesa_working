#   @(#)Makefile	1.4   8/2/91
#
# Makefile for TIMLIB
#
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
FC = fort
CPPFLAGS = -DDECPOINTER
FFLAGS = -c -O4 $(CPPFLAGS)
MAIN = /usr/home/bis/mesa/bmesa/Time_Propagation
TIMLIB = $(MAIN)/timlib.a

GET = sccs get

SRCS = \
	band3.F \
        band5.F \
	band7.F \
	chpsi.F \
	cpsi0.F \
	cpy3.F \
	cvect0.F \
	eigen.F \
	fdiff.F \
	gpaket.F \
	moment.F \
	nrpaket.F \
	packh.F \
	pert.F \
	points.F \
	scdip.F \
	setso.F \
	vcple.F \
	vnl.F \
	vtim.F \
	zfil2.F \
	zfil3.F \
	zproj.F

.F.o:
	$(FC) $(FFLAGS) $<

all: $(TIMLIB)
$(TIMLIB): $(SRCS:.F=.o)
	$(AR) $(ARFLAGS) $(TIMLIB) $(?:.F=.o)
