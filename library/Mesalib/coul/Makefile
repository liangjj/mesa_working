#    %W%   %G%
#
#    Makefile for the coulomb routine
#    
MAINDIR = /usr/home/bis/mesa
COULDIR = $(MAINDIR)/library/coul
MYLIB = $(MAINDIR)/mylib/mylib.a
AR = -ar
RANLIB = ranlib
FC = fort
FFLAGS =
REL = 

SRCS = \
	cf1a.f \
	cf1c.f \
	cf1r.f \
	cf2.f \
	clogam.f \
	coulcc.f \
	f11.f \
	f20.f \
	grncal.f \
	rcf.f

.f.o:
	$(F77) $(F77FLAGS) -c $<

all: $(MYLIB)
$(MYLIB): $(SRCS:.f=.o)
	$(AR) rvu $(MYLIB) $(?:.f=.o)
	$(RANLIB) $(MYLIB)


