include DEFS

PROG =	pgd

SRCS =	blacs.f90 error_prt.f90 pdsygvd2.f90 pg_diag.f90 \
	rd_oz_file.f90 precisn.f90 io_units.f90 pgd.f90 \
	nmlst.f90

OBJS =	$(SRCS:.f90=.o)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) $<

io_units.o: io_units.f90
precisn.o: precisn.f90
blacs.o: blacs.f90 io_units.o
error_prt.o: error_prt.f90 io_units.o
pdsygvd2.o: pdsygvd2.f90 blacs.o precisn.o
pg_diag.o: pg_diag.f90 blacs.o error_prt.o pdsygvd2.o rd_oz_file.o
rd_oz_file.o: rd_oz_file.f90 blacs.o nmlst.o io_units.o
pgd.o: pgd.f90 precisn.o io_units.o blacs.o pg_diag.o nmlst.o
nmlst.o: nmlst.f90 blacs.o io_units.o
