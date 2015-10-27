#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/CC_Prop
include $(TOP_DIR)/Makefile.inc
SRCS = \
	input_output.f90 \
	potential.f90 \
	grid_global.f90 \
	dvr_shared.f90 \
	dvr_global.f90 \
	prop_prnt.f90 \
	prop_global.f90 \
	Iterative_Global.f90\
	dvr_prnt.f90 \
	fd_prnt.f90 \
	fd_global.f90 \
	dvrprop_global.f90 \
	Channel_Module.f90


%.mod : %.o
	@if [! -f $@ ]; then \
	rm $< \
	$(MAKE) $< \
        fi

%.o : %.f90
	$(F90) $(F90FLAGS) -I$(MODULES) $<

all: $(MOD_LIB_A) 
$(MOD_LIB_A): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(MOD_LIB_A) $(?:.f90=.o)
sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod

