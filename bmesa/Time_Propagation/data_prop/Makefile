#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the arnoldi code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
MAIN = /usr/home/bis/mesa
LIB = $(MAIN)/library
MESABIN =  /usr/home/bis/mesa/bin
MESALIB = $(LIB)/mesalib.a

FC = fort
FFLAGS = -c -O4 

LD = $(FC)
LDFLAGS =

RCS = RCS
CO = co -l

SRCS = \
	data_prop.f90

.f90.o:
	$(FC) $(FFLAGS) -I$(MODDIR) $<

all: $(MESABIN)/data_prop
$(MESABIN)/data_prop: $(SRCS:.f90=.o)  \
                      $(MESALIB)

	$(LD) $(LDFLAGS) $(SRCS:.f90=.o)  \
	$(MESALIB) -o $(MESABIN)/data_prop     


sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm *.o

