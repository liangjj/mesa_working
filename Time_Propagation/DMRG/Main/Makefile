#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
#
RCS = RCS
CO = co -l
DMRG_MODULE_DIR = $(MESA_HOME)/Time_Propagation/DMRG/Modules
DMRG_MOD_LIB =$(DMRG_MODULE_DIR)/dmrg_mod_lib_f90.a
DMRG_DRIVER = $(MESA_HOME)/Time_Propagation/DMRG/Main
CPPFLAGS = MPI
CPPSRCS = \
	MPDFortran.F90 

.F90.o:
	$(F90) $(F90FLAGS) $(CPPFLAGS) -I$(DMRG_MODULE_DIR) $<

all: $(DMRG_DRIVER)/dmrg
$(DMRG_DRIVER)/dmrg: $(CPPSRCS:.F90=.o) \
			$(DMRG_MOD_LIB) $(MD_LIB)
			$(LD90) $(LDFLAGS) $(CPPSRCS:.F90=.o) \
			$(DMRG_MOD_LIB) $(MD_LIB) -o $(DMRG_DRIVER)/dmrg

files: $(SRCS)
	$(CO) $(RCS)/*.F90

clean:
	rm -f *.o 
