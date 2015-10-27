#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the potential library
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
RANLIB = ranlib
REL = 

SRCS = \
	addang.f90 \
	mwell.f90 \
	none.f90 \
	rwell.f90 \
	vadex12.f90 \
	vadpr.f90 \
	vanhmo.f90 \
	vcad12.f90 \
	vcoul.f90 \
	v_eberlonium.f90 \
	vexp.f90 \
	vexp_sum.f90 \
	vexp12.f90 \
	vexp123.f90 \
	vhmo.f90 \
	vir4.f90 \
	vmorse.f90 \
	vone.f90 \
	vperiod.f90 \
	vpert.f90 \
	vres.f90 \
	vrwell.f90 \
	v_temkin_poet.f90 \
	vwad12.f90 \
	vwell.f90 \
	vwel12.f90 \
	vwel123.f90 \
	vyukawa.f90 \
	v_pow_exp.f90 \
	v_e_laser.f90

.f90.o:
	$(F90) $(F90FLAGS) -I$(MOD_DIR) $<

all: $(POTENTIAL_LIB)
$(POTENTIAL_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(POTENTIAL_LIB) $(?:.f90=.o)

sources: $(SRCS)
	co -l RCS/*.f90
clean:
	rm -f *.o
	rm -f *.a
