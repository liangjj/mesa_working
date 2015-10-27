#@(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR = $(MESA_HOME)/Legendre
include $(TOP_DIR)/Makefile.inc
SRCS = \
	Driver.f90 

all: $(DRIVER) 
MD_LIB=$(MD_LIB_SEQUENTIAL_DYNAMIC)
$(DRIVER): $(SRCS:.f90=.o) $(LIB)
	$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) $(LIB) $(MD_LIB) \
	$(GOTOBLASLIB) -o $(DRIVER)



files: $(SRCS) 
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
