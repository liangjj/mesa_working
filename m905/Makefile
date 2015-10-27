#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the arnoldi code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
MAIN = /usr/home/bis/mesa
LIB = $(MAIN)/library
MESABIN =  /usr/home/bis/mesa/bin
MESALIB = $(LIB)/mesalib.a
MDLIB = /usr/lib/libcxml.a    
MODDIR = $(LIB)/Modules
MODLIB = $(LIB)/modlib_f90.a

FC = fort
FFLAGS = -c -O4 

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	drive_lanczos.f90 \
	init_vec.f90 \
	lanczos.f90 \
	lanczos_mult.f90 \
	rdmat.f90

.f90.o:
	$(FC) $(FFLAGS) -I$(MODDIR) $<

all: $(MESABIN)/m905
$(MESABIN)/m905: $(SRCS:.f90=.o)  \
                 $(MODLIB) $(MESALIB) $(MDLIB)
		 $(LD) $(LDFLAGS) $(SRCS:.f90=.o)  \
	         $(MODLIB) $(MESALIB) $(MDLIB) -o $(MESABIN)/m905     


sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o

