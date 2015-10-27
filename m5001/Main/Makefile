# @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
SOURCE = m5001
include $(MESA_HOME)/$(SOURCE)/Makefile.inc
SRCLOC = $(MESA_HOME)/$(SOURCE)/Main
SRCS =	$(SRCLOC)/$(SOURCE).f90
all: $(MESA_BIN)/$(SOURCE)
$(MESA_BIN)/$(SOURCE): $(SRCS:.f90=.o) 
	$(F90) \
	$(SRCS:.f90=.o)	\
	$(LOAD_LIB) $(MD_LIB) -o $(MESA_BIN)/$(SOURCE)  

clean:
	rm -f *.o *.mod *~


#	$(LOAD_F90) -o $(MESA_BIN)/$(SOURCE)  