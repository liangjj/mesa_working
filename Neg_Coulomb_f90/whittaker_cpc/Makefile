#!  Time-stamp: "2001-05-15 15:03:31 cjn"
PROG1 =	whittaker
PROG2 = simos

SRCS =	precisn.f90 \
	io_units.f90 \
	weniger.f90 \
	whittaker_w.f90

SRCS1 = whittaker_w_test.f90
SRCS2 = simos.f90 simos_test.f90

OBJS =	precisn.o \
	io_units.o \
	weniger.o \
	whittaker_w.o

OBJ1 = whittaker_w_test.o
OBJ2 = simos.o simos_test.o

#IBM RS6000:
#LIBS = -L/usr/tci26/lib -lessl_r -llapack
#F90 = xlf95_r
#F90FLAGS = -F:f90 -g
#F90FLAGS = -F:f90 -O5 -qarch=pwr3 -qtune=pwr3 -qcache=auto -qnosave \
	-qlanglvl=95pure -qthreaded -qreport=hotlist -qipa=level=2
#LDFLAGS =

#SUN Ultra-5_10 SPARC; Solaris (SunOS 5.8)
#LIBS = -xlic_lib=sunperf
#F90 = f95
#F90FLAGS = -c -dalign -fast -depend
#LDFLAGS = -dalign -fast -depend

#SGI Origin2000
#LIBS = -lcomplib.sgimath
#F90 = f90
#F90FLAGS = -Ofast=ip27 -64 -OPT:Olimit=0
#LDFLAGS = -64 -IPA

#Cray T3E-1200E
#F90 = f90
#F90FLAGS = -O3 -e n
#LDFLAGS = -Wl"-Dequiv(DSTEQR)=SSTEQR"
#LIBS = 

#Compaq alpha; Redhat Linux 6.1
#F90 = fort
#F90FLAGS = -O -fast -tune host -arch host -check overflow -std95
#LDFLAGS =
#LIBS = -lcxml 


all: $(PROG1) $(PROG2)

$(PROG1): $(OBJS) $(OBJ1)
	$(F90) $(LDFLAGS) $(OBJS) $(OBJ1) $(MD_LIB) -o $@
$(PROG2): $(OBJS) $(OBJ2)
	$(F90) $(LDFLAGS) $(OBJ2) $(OBJS) $(MD_LIB) -o $@

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

weniger.o: precisn.o io_units.o
whittaker_w.o: precisn.o io_units.o weniger.o
simos_test.o: precisn.o io_units.o whittaker_w.o simos.o
whittaker_w_test.o: precisn.o io_units.o whittaker_w.o











