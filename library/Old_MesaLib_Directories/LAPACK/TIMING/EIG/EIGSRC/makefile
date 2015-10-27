FORTRAN = f77
OPTS    = -O -u
NOOPT   = -u

#######################################################################
#  This is the makefile to create a library of the eigensystem routines
#  from LAPACK that have been instrumented to count operations. 
#  The files are organized as follows:
#
#     SCIGSRC -- Single precision instrumented LAPACK routines used in
#                REAL and COMPLEX
#     DZIGSRC -- Double precision instrumented LAPACK routines used in
#                DOUBLE PRECISION and COMPLEX*16
#
#     SEIGSRC -- Single precision real instrumented LAPACK routines
#     CEIGSRC -- Single precision complex instrumented LAPACK routines
#     DEIGSRC -- Double precision real instrumented LAPACK routines
#     ZEIGSRC -- Double precision complex instrumented LAPACK routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  First, modify the FORTRAN and OPTS 
#  definitions to match your compiler and the options to be used.
#  Then to create or add to the library, enter make followed by one or
#  more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       eigsrc.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#----------------------------------------------------------------------
#
#  Edward Anderson, University of Tennessee
#  March 26, 1990.
#  Susan Ostrouchov, Last updated January 19, 1993.
#
#######################################################################

SCIGSRC = slaebz.o sstebz.o

DZIGSRC = dlaebz.o dstebz.o
 
SEIGSRC = \
   sbdsqr.o sgghrd.o shgeqz.o shsein.o shseqr.o \
   slaein.o \
   slahqr.o spteqr.o sstein.o ssteqr.o ssterf.o \
   stgevc.o strevc.o

CEIGSRC = \
   cbdsqr.o cgghrd.o chgeqz.o chsein.o chseqr.o \
   claein.o clahqr.o \
   cpteqr.o cstein.o csteqr.o ctgevc.o ctrevc.o

DEIGSRC = \
   dbdsqr.o dgghrd.o dhgeqz.o dhsein.o dhseqr.o \
   dlaein.o \
   dlahqr.o dpteqr.o dstein.o dsteqr.o dsterf.o \
   dtgevc.o dtrevc.o

ZEIGSRC = \
   zbdsqr.o zgghrd.o zhgeqz.o zhsein.o zhseqr.o \
   zlaein.o zlahqr.o \
   zpteqr.o zstein.o zsteqr.o ztgevc.o ztrevc.o

all: single complex double complex16

single: $(SEIGSRC) $(SCIGSRC)
	ar cr ../eigsrc.a $(SEIGSRC) $(SCIGSRC)
	ranlib ../eigsrc.a

complex: $(CEIGSRC) $(SCIGSRC)
	ar cr ../eigsrc.a $(CEIGSRC) $(SCIGSRC)
	ranlib ../eigsrc.a

double: $(DEIGSRC) $(DZIGSRC)
	ar cr ../eigsrc.a $(DEIGSRC) $(DZIGSRC)
	ranlib ../eigsrc.a

complex16: $(ZEIGSRC) $(DZIGSRC)
	ar cr ../eigsrc.a $(ZEIGSRC) $(DZIGSRC)
	ranlib ../eigsrc.a

$(SCIGSRC): $(FRC)
$(DZIGSRC): $(FRC)
$(SEIGSRC): $(FRC)
$(CEIGSRC): $(FRC)
$(DEIGSRC): $(FRC)
$(ZEIGSRC): $(FRC)

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o

.f.o: 
	$(FORTRAN) -c $(OPTS) $<
