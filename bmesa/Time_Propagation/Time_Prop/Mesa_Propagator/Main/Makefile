#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
MESA_PROPAGATOR_DIR = ..
MESA_LIB_DIR =$(MESA_LIB)
MESA_LIBRARY =$(MESA_LIB_DIR)/mesalib.a
MOD_DIR = $(MESA_LIB_DIR)/Modules
MOD_LIB = $(MOD_DIR)/modlib_f90.a
#DVR_DIR = $(MESA_LIB_DIR)/dvrlib_f90
#DVR_LIB = $(DVR_DIR)/dvrlib_f90.a
#FD_DIR = $(MESA_LIB_DIR)/fdlib_f90
#FD_LIB = $(FD_DIR)/fdlib_f90.a
DVR_DIR = $(MESA_LIB)/DVR_Library
DVR_LIB = $(DVR_DIR)/dvr_lib_f90.a
FD_LIB = $(DVR_DIR)/fd_lib_f90.a
POT_DIR = $(MESA_LIB_DIR)/potlib
POT_LIB = $(POT_DIR)/potlib_f90.a
PROP_MOD_DIR = $(MESA_PROPAGATOR_DIR)/Prop_Modules
PROP_SOURCES_DIR = $(MESA_PROPAGATOR_DIR)/Prop_Sources
PROP_SOURCE_LIB = $(PROP_SOURCES_DIR)/prop_source_lib_f90.a
PROP_MOD_LIB = $(PROP_MOD_DIR)/prop_mod_lib_f90.a

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

%.mod : %.o
	@if [! -f $@ ]; then \
	rm $< \
	$(MAKE) $< \
	fi

%.o : %.f90
	$(FC) $(FFLAGS) -I$(MOD_DIR) -I$(PROP_MOD_DIR) $<

all: $(MESA_BIN)/m1235
$(MESA_BIN)/m1235: $(SRCS:.f90=.o) \
                  $(PROP_MOD_LIB) $(PROP_SOURCE_LIB) \
                  $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POT_LIB) \
                  $(MESA_LIBRARY) $(MD_LIB)
	$(LD) $(LDFLAGS) $(SRCS:.f90=.o) \
        $(PROP_MOD_LIB) $(PROP_SOURCE_LIB) \
        $(DVR_LIB) $(FD_LIB) $(MOD_LIB) $(POT_LIB)  \
	$(MESA_LIBRARY) $(MD_LIB) -o $(MESA_BIN)/m1235     



files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o 
