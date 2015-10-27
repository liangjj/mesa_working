#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
#
DVR_DRIVER = $(B_MESA)/DVR_Driver
DVR_DRIVER_MAIN = $(DVR_DRIVER)/DVR_Main
DVR_DRIVER_MOD = $(DVR_DRIVER)/DVR_Modules
DVR_DRIVER_SOURCES = $(DVR_DRIVER)/DVR_Sources
#
DVR_DRIVER_LIB = $(DVR_DRIVER)/dvr_driver_lib_f90.a
#

RCS = RCS
CO = co -l

SRCS = dvr_driver.f90

F90FLAGS = -L$(MKL_PATH) -I$(MKL_INCLUDE) -I$(MODULES) -I$(DVR_DRIVER_MOD) -O4
LDFLAGS = $(F90FLAGS)

.f90.o:
	$(F90) $(F90FLAGS) -c $<

all: $(DVR_DRIVER_MAIN)/dvr_x
$(DVR_DRIVER_MAIN)/dvr_x: $(SRCS:.f90=.o) \
			$(DVR_DRIVER_LIB) \
			$(DVR_LIB) \
			$(FD_LIB) \
			$(MOD_LIB) \
			$(POTENTIAL_LIB) \
			$(UTILITY_LIB) 
	$(LD90) $(LDFLAGS) $(SRCS:.f90=.o) \
			$(DVR_DRIVER_LIB) \
			$(DVR_LIB) \
			$(FD_LIB) \
			$(MOD_LIB) \
			$(POTENTIAL_LIB)  \
			$(UTILITY_LIB) \
			$(MD_LIB) \
	-o $(DVR_DRIVER_MAIN)/dvr_x     



files: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o 
