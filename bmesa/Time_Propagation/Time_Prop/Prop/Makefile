#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
REMOVE = NO
MAIN = $(MESA_HOME)
LIB = $(MAIN)/library
LIBMODDIR = $(LIB)/Modules
MODDIR = Modules
SOURCEDIR = Sources
MESABIN =  $(MESA_BIN)
MODLIB = $(LIB)/modlib_f90.a
SOURCELIB = Sources/sourcelib.a
PROPMODLIB = Modules/propmodlib.a
TIMPRP = $(MAIN)/bmesa/Time_Propagation
THISDIR = $(TIMPRP)/Time_Prop/Prop
DVRLIB = $(LIB)/dvrlib_f90/dvrlib_f90.a
FDLIB = $(LIB)/fdlib_f90/fdlib_f90.a
POTLIB = $(LIB)/potlib/potlib_f90.a
MESALIB = $(LIB)/mesalib.a
ITDVRLIB  = $(LIB)/itlib/itdvrlib.a
TDVRLIB = $(LIB)/tdvrlib/tdvrlib.a

FFLAGS = -c -O5 -I$(LIBMODDIR) -I$(MODDIR)
LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = prop_main.f90

.f90.o:
	$(FC) $(FFLAGS) $<

all: $(MESABIN)/m1235
$(MESABIN)/m1235: $(SRCS:.f90=.o) \
                  $(PROPMODLIB) $(SOURCELIB) \
                  $(DVRLIB) $(FDLIB) \
                  $(MODLIB) $(POTLIB) $(ITDVRLIB) \
                  $(TDVRLIB) $(MESALIB) $(MDLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f90=.o) \
        $(PROPMODLIB) $(SOURCELIB) \
	$(DVRLIB) $(FDLIB) $(MODLIB) \
	$(POTLIB) $(ITDVRLIB) $(TDVRLIB) \
	$(MESALIB) $(MDLIB) -o $(MESABIN)/m1235     



files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o *.mod *~
