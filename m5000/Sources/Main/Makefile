#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/m5000
include $(MESA_HOME)/Makefile.inc
include $(TOP_DIR)/Makefile.inc
SRCS =	\
	m5000.f90
#
LOAD =  $(LOAD_f90) $(LOAD_LIB)
all: $(MESA_BIN) $(SRCS:.f90=.o) 
	$(LOAD) -o $(MESA_BIN)/m5000     

clean:
	rm -f *.o *.mod *~
