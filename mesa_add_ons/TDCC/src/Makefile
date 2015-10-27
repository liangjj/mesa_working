FC = mpif90
LD = $(FC)

F90C = $(FC)
F90FLAGS = $(FFLAGS)
CXXFLAGS = $(FFLAGS)

.SUFFIXES: .f90 .mod $(SUFFIXES)

DEFAULT: TDCC_IT.exe

PHONY: clean

# .f files are Fortran77 and do not contain modules
%.o: %.f
	$(FC) $(FFLAGS) -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $<

%.mod: %.f90
	$(F90C) $(F90FLAGS) -c $<

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $<

%.mod: %.F90
	$(F90C) $(F90FLAGS) -c $<

%.exe : %.o
	$(F90C) $(F90FLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm -f *.o *.mod *.il *.exe

### dependencies as found by makeMakefile.py
anglib.o anglib.mod : nrtype.mod

invmod.o invmod.mod : nrtype.mod

pfedvrmod.o pfedvrmod.mod : nrtype.mod

singleel.o singleel.mod : dvrmod.mod nrtype.mod pfedvrmod.mod

TDCC_IT.o : anglib.mod dvrmod.mod globalmod.mod integrals.mod invmod.mod keops.mod nrtype.mod pfedvrmod.mod tdcc_itsubroutines.mod

TDCC_IT.exe : anglib.o dvrmod.o globalmod.o integrals.o invmod.o keops.o nrtype.o pfedvrmod.o tdcc_itsubroutines.o

# slagge.o zlagge.o dlagge.o clagge.o

singleeltest.o : blas95.mod dvrmod.mod nrtype.mod pfedvrmod.mod singleel.mod

singleeltest.exe : blas95.o dvrmod.o nrtype.o pfedvrmod.o singleel.o

keops.o keops.mod : globalmod.mod pfedvrmod.mod

globalmod.o globalmod.mod : nrtype.mod

tdcc_itsubroutines.o tdcc_itsubroutines.mod : dvrmod.mod globalmod.mod pfedvrmod.mod

