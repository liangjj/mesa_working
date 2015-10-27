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
TIMPRP = $(MAIN)/bmesa/Time_Propagation
EIGENPROP = $(TIMPRP)/Eigen_Prop_f90
DVRLIB = $(LIB)/dvrlib_f90/dvrlib_f90.a
FDLIB = $(LIB)/fdlib_f90/fdlib_f90.a
POTLIB = $(LIB)/potlib/potlib_f90.a
HAMLIB = $(LIB)/hamlib/hamlib.a
MESALIB = $(LIB)/mesalib.a
ITDVRLIB  = $(LIB)/itlib/itdvrlib.a
TDVRLIB = $(LIB)/tdvrlib/tdvrlib.a
TIMLIB = $(TIMPRP)/timlib_f90.a
TPRPLIB = $(TIMPRP)/tprplib_f90.a
MDLIB = /usr/lib/libcxml.a    

FC = fort
FFLAGS = -c -O4 

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	eigen_main.f90 \
	eigen_prop.f90 

.f90.o:
	$(FC) $(FFLAGS) -I$(MODDIR) $<

all: $(MESABIN)/m1110
$(MESABIN)/m1110: $(SRCS:.f90=.o)  \
                  $(TIMLIB) $(TPRPLIB) $(DVRLIB) $(FDLIB) \
                  $(MODLIB) $(HAMLIB) $(POTLIB) $(ITDVRLIB) \
                  $(TDVRLIB) $(MESALIB) $(MDLIB)

	$(LD) $(LDFLAGS) $(SRCS:.f90=.o)  \
	$(TIMLIB) $(TPRPLIB) $(DVRLIB) $(FDLIB) $(MODLIB) \
	$(HAMLIB) $(POTLIB) $(ITDVRLIB) $(TDVRLIB) \
	$(MESALIB) $(MDLIB) -o $(MESABIN)/m1110     


sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o

