#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
TOP_DIR = $(MESA_LIB)/Spheroidal_DVR
include $(TOP_DIR)/Makefile.inc
RCS = RCS
CO = co -l
REL = 

SRCS = \
	FEDVR_Hamiltonian_Module.f90

all: $(SPHEROIDAL_DVR_LIB)
$(SPHEROIDAL_DVR_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(SPHEROIDAL_DVR_LIB) $(?:.f90=.o)
	$(RANLIB) $(SPHEROIDAL_DVR_LIB)
sources: $(SRCS)
	$(CO) $(RCS)/*.f90*

clean:
	rm -f *.o

