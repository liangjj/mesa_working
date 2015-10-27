#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the arnoldi code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
MAIN = /usr/home/bis/mesa
LIB = $(MAIN)/library
MODDIR = $(LIB)/Modules
MESABIN =  /usr/home/bis/mesa/bin
MODLIB = $(LIB)/modlib_f90.a
DVRLIB = $(LIB)/dvrlib_f90/dvrlib_f90.a
POTLIB = $(LIB)/potlib/potlib_f90.a
MESALIB = $(LIB)/mesalib.a
MDLIB = /usr/lib/libcxml.a    

FC = fort
FFLAGS = -c -O4 

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	r_matrix_main.f90 \
	r_mat.f90 
.f90.o:
	$(FC) $(FFLAGS) -I$(MODDIR) $<

all: $(MESABIN)/m1102
$(MESABIN)/m1102: $(SRCS:.f90=.o)  \
                    $(DVRLIB) $(MODLIB) $(POTLIB)  \
                    $(MESALIB) $(MDLIB)

	$(LD) $(LDFLAGS) $(SRCS:.f90=.o)  \
	$(DVRLIB) $(MODLIB) $(POTLIB)  \
	$(MESALIB) $(MDLIB) -o $(MESABIN)/m1102     


sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o

