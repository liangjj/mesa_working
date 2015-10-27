include DEFS
PROG =	DC_DIAG
MODSDIR = 
BIN_DIR = $(MESA_HOME)/Test_DC/BIN_DIR
WORK_DIR = $(MESA_HOME)/Test_DC/WORK_DIR
ZCOM = Read_ipar.f90 Read_apar.f90 Read_rpar.f90 Check_file.f90 

SRCS =	mod_io_units.f90 mod_precisn.f90 mod_blacs.f90  mod_dc_matrix.f90  \
	pgd_main.f90 pgd_drive.f90   pdsygvd1.f90 error_prt.f90 read_dcmat.f90 write_evec.f90

OBJS =	$(SRCS:.f90=.o) $(ZCOM:.f90=.o)
F90=mpif90
F90FLAGS=-O3
LDFLAGS=-L/opt/Intel/composer_xe_2011_sp1.7.256/mkl/lib/intel64 -I/opt/Intel/composer_xe_2011_sp1.7.256/mkl/include
MD_LIB_MPI=-L/opt/Intel/composer_xe_2011_sp1.7.256/mkl/lib/intel64 -lmkl_scalapack_lp64 /opt/Intel/composer_xe_2011_sp1.7.256/mkl/lib/intel64/libmkl_solver_lp64.a -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -Wl,--end-group -openmp -lpthread -lm
#
#    When using openmpi you need to link statically using
NOLINK = -c
all: $(PROG)  

$(PROG): $(OBJS)
#	$(F90) $(F90FLAGS)  $(LDFLAGS) $(OBJS) $(UTILITY_LIB) $(MD_LIB_MPI) -o $@
	$(F90) $(F90FLAGS)  $(LDFLAGS) $(OBJS) $(MD_LIB_MPI) -o $@
	cp $(PROG) $(BIN_DIR)
	cp $(PROG) $(WORK_DIR)
	rm *.o *.mod

clean:
	rm -f $(OBJS)  *.mod *~

x:
	rm -f $(PROG)

.SUFFIXES:
.SUFFIXES: .f90 .f77 .mod .o

.f90.o:
	$(F90) $(NOLINK) $(F90FLAGS) $<
.f.o:
	$(F77) $(NOLINK) $(FFLAGS) $<

