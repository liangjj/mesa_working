#@(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR = $(MESA_HOME)/Test_Array
include $(TOP_DIR)/Makefile.inc
SRCS = \
	Driver.f90 

all: $(DRIVER) 
$(DRIVER): $(SRCS:.f90=.o) $(LIB)
	$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) $(LOAD_LIB) $(MD_LIB) \
	-o $(DRIVER)



files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
