
#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the dvr code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
MAIN = /usr/home/bis/mesa
LIB = $(MAIN)/library
MODDIR = $(LIB)/Modules
MESABIN =  /usr/home/bis/mesa/bin
MODLIB = $(LIB)/modlib_f90.a
DVR = $(MAIN)/bmesa/Time_Independent
DVRLIB = $(LIB)/dvrlib_f90/dvrlib_f90.a
FDLIB = $(LIB)/fdlib_f90/fdlib_f90.a
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
	dvr_main.f90 \
	dvd_dat.f90 \
	v_couple.f90

.f90.o:
	$(FC) $(FFLAGS) -I$(MODDIR) $<

all: $(MESABIN)/m1109
$(MESABIN)/m1109: $(SRCS:.f90=.o)  \
                    $(DVRLIB) $(FDLIB) \
                    $(MODLIB) $(POTLIB) \
                    $(MESALIB) $(MDLIB)

	$(LD) $(LDFLAGS) $(SRCS:.f90=.o)  \
	$(DVRLIB) $(FDLIB) $(MODLIB) \
	$(POTLIB)  \
	$(MESALIB) $(MDLIB) -o $(MESABIN)/m1109     


sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o

