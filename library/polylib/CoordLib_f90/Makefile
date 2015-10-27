#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
COORDLIB = ../coordlib_f90.a
MODDIR =/usr/home/bis/mesa/library/Modules
.SUFFIXES: .f90 $(SUFFIXES)
AR = ar
ARFLAGS = rvu
FC = fort

FFLAGS = -c -O4 -I$(MODDIR)
LD =$(FC)
LDFLAGS =
RCS = RCS
CO = co -l
REL = 

SRCS = \
	cordfn.f90 \
	cpoly.f90 \
	drvply.f90 \
	lgngr.f90

.f90.o:
	$(FC) $(FFLAGS) $<

all: $(COORDLIB)
$(COORDLIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(COORDLIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90*

clean:
	rm *.o

