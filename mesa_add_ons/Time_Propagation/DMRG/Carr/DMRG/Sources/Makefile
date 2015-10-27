#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = ..
include $(TOP_DIR)/Makefile.inc

MOD_INCLUDE = -I$(DMRG_MODULE_DIR)

SRCS = \
	cebct.f90 \
	ebct.f90 

all: $(DMRG_SOURCES_LIB)

$(DMRG_SOURCES_LIB): $(SRCS:.f90=.o) 
	$(AR) $(ARFLAGS) $(DMRG_SOURCES_LIB) $(?:.f90.o)
	$(RANLIB) $(DMRG_SOURCES_LIB)

files: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o 
