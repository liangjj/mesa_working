#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
GENERIC_PROPAGATOR_DIR = \
/home/mesa/mesa/bmesa/Time_Propagation/Time_Prop/Generic_Propagator
MOD_DIR = $(GENERIC_PROPAGATOR_DIR)/Modules
MOD_LIB = $(MOD_DIR)/mod_lib_f90.a
DVR_LIBRARY =  \
$(GENERIC_PROPAGATOR_DIR)/DVR_Library
DVR_LIB = $(DVR_LIBRARY)/dvr_lib_f90.a
FD_LIB = $(DVR_LIBRARY)/fd_lib_f90.a
PROP_MOD_DIR = Prop_Modules
PROP_SOURCES_DIR = Prop_Sources
PROP_SOURCE_LIB = $(PROP_SOURCES_DIR)/prop_source_lib_f90.a
PROP_MOD_LIB = $(PROP_MOD_DIR)/prop_mod_lib_f90.a
THISDIR = ../Prop_To_Go
BIN = $(THISDIR)
UTIL_LIB = ../Utility_To_Go/util_lib_f90.a
POTLIB = ../Potential/pot_lib_f90.a
LIBPATH = /usr/pgi/linux86-64/6.1/lib

#LIBMKL_EM64T = /opt/intel/mkl721/lib/em64t/libmkl_em64t.a 
#LIBMKL_LAPACK = /opt/intel/mkl721/lib/em64t/libmkl_lapack.a 
#LIBGUIDE = /opt/intel/mkl721/lib/em64t/libguide.a 
#LIBMKL_SOLVER = /opt/intel/mkl721/lib/em64t/libmkl_solver.a 
#LIBPTHREAD=/usr/lib64/libpthread.a 
#LIBC=/usr/lib64/libc.a
#MDLIB = $(LIBMKL_LAPACK)$(LIBMKL_SOLVER) $(LIBMKL_EM64T) $(LIBGUIDE) $(LIBPTHREAD) \
        $(LIBC)
LD = $(FC)
LDFLAGS = 

RCS = RCS
CO = co -l

SRCS = prop_main.f90

.f90.o:
	$(FC) $(FFLAGS) -I$(MOD_DIR) -I$(PROP_MOD_DIR) $<

all: $(BIN)/dvrprop
$(BIN)/dvrprop: $(SRCS:.f90=.o) \
                  $(PROP_MOD_LIB) $(PROP_SOURCE_LIB) \
                  $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POTLIB) \
                  $(UTIL_LIB) $(MD_LIB)
	$(LD) $(LDFLAGS) $(SRCS:.f90=.o) \
        $(PROP_MOD_LIB) $(PROP_SOURCE_LIB) \
        $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POTLIB)  \
	$(UTIL_LIB) $(MD_LIB) -o $(BIN)/dvrprop     



files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o 
