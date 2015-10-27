# Time-stamp: "2008-11-02 04:27:58 Clifford Noble"
#include DEFS_tacc
include DEFS_intel
#include DEFS_bis

SRCS =	blacs.f90 gen_diag.f90 pg_diag.f90 \
	rd_dcmat.f90 precisn.f90 io_units.f90 pgd.f90 \
	nmlst.f90

OBJS =	$(SRCS:.f90=.o)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod *.o

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c  $<

io_units.o: io_units.f90
precisn.o: precisn.f90
blacs.o: blacs.f90 io_units.o
gen_diag.o: gen_diag.f90 blacs.o precisn.o
pg_diag.o: pg_diag.f90 blacs.o gen_diag.o rd_dcmat.o
rd_dcmat.o: rd_dcmat.f90 blacs.o nmlst.o io_units.o
pgd.o: pgd.f90 precisn.o io_units.o blacs.o pg_diag.o nmlst.o
nmlst.o: nmlst.f90 blacs.o io_units.o
