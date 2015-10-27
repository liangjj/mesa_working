#   @(#)Makefile	1.4   8/2/91
#

TOP_DIR = $(MESA_HOME)/Intel_Test
#
# Suffix rules
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
#
.F90.o:
	$(F90) $(F90FLAGS) $<
.f90.o:
	$(F90) $(F90FLAGS) $<

SRCS =	barry.f90
#
DRIVER = $(TOP_DIR)/cliff
#
all: $(DRIVER)
$(DRIVER): $(SRCS:.f90=.o) $(MD_LIB)
	$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) $(MD_LIB) \
	-o $(DRIVER)     
