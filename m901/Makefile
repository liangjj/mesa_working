#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M901
#
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
FC = $(FC)
FFLAGS = $(FFLAGS)

LD = fort
LDFLAGS =

BINDIR = $(MESA_BIN)
MESALIB = $(MESA_LIB)/mesalib.a
MDLIB = $(MDLIB)
M901LIB = $(MESA_LIB)/m901lib.a
GET = 
REL = 

SRCS = \
	../m901/m901.f \
	../m901/loopdat.f \
	../m901/pm901.f

LIBSRCS = \
	../m901/fdiag.f \
	../m901/fold.f \
	../m901/foldw.f \
	../m901/getdrt.f \
	../m901/inhmat.f \
	../m901/init4x.f \
	../m901/initdg.f \
	../m901/initex.f \
	../m901/initlp.f \
	../m901/initpg.f \
	../m901/initst.f \
	../m901/initx2.f \
	../m901/loopdat.f \
	../m901/loopyd.f \
	../m901/mn901.f \
	../m901/roots.f \
	../m901/square.f \
	../m901/squarw.f \
	../m901/start.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m901
lib: $(M901LIB)
$(BINDIR)/m901: $(SRCS:.f=.o) $(M901LIB) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(M901LIB) $(MESALIB) \
                $(MDLIB) -o $(BINDIR)/m901
$(M901LIB): $(LIBSRCS:.f=.o)
	$(AR) $(ARFLAGS) $(M901LIB) $(?:.f=.o)
	$(RANLIB) $(M901LIB)

clean:
	rm -f *.o core
