#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
GENERIC_PROPAGATOR_DIR = \
$(MESA_HOME)/bmesa/Time_Propagation/Time_Prop/Generic_Propagator
MOD_DIR = $(GENERIC_PROPAGATOR_DIR)/Modules
MOD_LIB = $(MOD_DIR)/mod_lib_f90.a
DVR_DRIVER_DIR = $(GENERIC_PROPAGATOR_DIR)/DVR_Driver
DVR_DRIVER_MAIN_DIR = $(DVR_DRIVER_DIR)/DVR_Main
DVR_DRIVER_MOD_DIR = $(DVR_DRIVER_DIR)/DVR_Modules
DVR_DRIVER_SOURCES_DIR = $(DVR_DRIVER_DIR)/DVR_Sources
DVR_DRIVER_MOD_LIB = $(DVR_DRIVER_MOD_DIR)/dvr_driver_modules_lib_f90.a
DVR_DRIVER_SOURCES_LIB = $(DVR_DRIVER_SOURCES_DIR)/dvr_driver_sources_lib_f90.a
DVR_LIBRARY =  \
$(GENERIC_PROPAGATOR_DIR)/DVR_Library
DVR_LIB = $(DVR_LIBRARY)/dvr_lib_f90.a
FD_LIB = $(DVR_LIBRARY)/fd_lib_f90.a
UTILITY_TO_GO_DIR = $(GENERIC_PROPAGATOR_DIR)/Utility_To_Go
UTIL_LIB = $(UTILITY_TO_GO_DIR)/util_lib_f90.a
POT_DIR = $(GENERIC_PROPAGATOR_DIR)/Potential
POT_LIB = $(POT_DIR)/pot_lib_f90.a
BIN = $(DVR_DRIVER_MAIN_DIR)
LIBPATH = /usr/pgi/linux86-64/6.1/lib
LD = $(FC)
LDFLAGS = 

RCS = RCS
CO = co -l

SRCS = dvr_driver.f90

.f90.o:
	$(FC) $(FFLAGS) -I$(MOD_DIR) -I$(DVR_DRIVER_MOD_DIR) $<

all: $(BIN)/dvr
$(BIN)/dvr: $(SRCS:.f90=.o) \
                  $(DVR_DRIVER_MOD_LIB) $(DVR_DRIVER_SOURCES_LIB) \
                  $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POT_LIB) \
                  $(UTIL_LIB) $(MD_LIB)
	$(LD) $(LDFLAGS) $(SRCS:.f90=.o) \
        $(DVR_DRIVER_MOD_LIB) $(DVR_DRIVER_SOURCES_LIB) \
        $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POT_LIB)  \
	$(UTIL_LIB) $(MD_LIB) -o $(BIN)/dvr     



files: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o 
