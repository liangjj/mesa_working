#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the mesa data code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
MAIN = /usr/home/bis/mesa
LIB = $(MAIN)/library
MODDIR = $(LIB)/Modules
MODLIB = $(LIB)/modlib_f90.a
MESABIN = $(MAIN)/bin
MESA_DATA_TRANSFER = $(MAIN)/m8000
MESALIB = $(LIB)/mesalib.a
MDLIB = /usr/lib/libcxml.a    

FC = fort
FFLAGS = -c -O5 -I$(MODDIR)

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	basout.f90 \
	mesa_data_transfer.f90 \
	newmat.f90
.f90.o:
	$(FC) $(FFLAGS) $<

all: $(MESABIN)/m8000
$(MESABIN)/m8000: $(SRCS:.f90=.o)  \
                  $(MESALIB) $(MODLIB) $(MDLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f90=.o)  \
	$(MESALIB) $(MODLIB) $(MDLIB) -o $(MESABIN)/m8000     

sources: $(SRCS) $(LIBSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o *.mod *~

#	mesa_data_transfer.f90 \
#	newmat.f90 \
#	trimc.f90