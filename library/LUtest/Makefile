#
#  Makefile for SuperLU demonstration
#  First Makes a version of SLUBridgelib.a in the local directory
#  Then compiles the fortran routines and links 
#
#  Assumes SuperLU standard installation in /usr/local/SuperLU
#
#  Linux version uses f77 and gcc
#

include /usr/local/SuperLU/make.inc

HEADER   = /usr/local/SuperLU/SRC

SOURCE= SLU_fortran.c
OBJECT=$(SOURCE:.c=.o)
BRIDGE=SLUBridgelib.a

DEFAULT: $(BRIDGE) $(BINDIR)/xtest

$(BRIDGE): $(OBJECT)
	ar cr $(BRIDGE) $(OBJECT)

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)


FC = f77
FFLAGS =  -c 

LD = f77
LDFLAGS =


SRCS = \
	  main.f matrix.f matxvec.f 

.f.o:
	$(FC)  -O2 -c  $<

$(BINDIR)/xtest: $(SRCS:.f=.o) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) -o $(BINDIR)xtest \
	./SLUBridgelib.a \
	/usr/local/SuperLU/superlu_linux.a \
	/usr/local/SuperLU/blas_linux.a \
	/usr/local/LAPACK/liblapack.a \
	/usr/local/LAPACK/libblas.a 

clean:
	rm -f *.o *.a xtest

