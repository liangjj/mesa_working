#   %W%   %G%
#
# Makefile for M820
#
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
#FC = g77
#FFLAGS = -c -g -finit-local-zero -g -finit-local-zero

LD = $(FC)
#LDFLAGS =

BINDIR = $(MESA_BIN)
MESALIB = $(MESA_LIB)/mesalib.a
M820LIB = $(MESA_LIB)/m820lib.a
GET = sccs get
REL = 

SRCS = \
	../m820/m820.f \
	../m820/pm820.f

LIBSRCS = \
	../m820/guglab.f \
	../m820/mn820.f \
	../m820/sriiii.f \
	../m820/sriijj.f \
	../m820/srijij.f \
	../m820/srijkl.f \
	../m820/srone.f \
	../m820/toguga.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m820
lib: $(M820LIB)
$(BINDIR)/m820: $(SRCS:.f=.o) $(M820LIB) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(M820LIB) $(MESALIB)\
                 -o $(BINDIR)/m820
$(M820LIB): $(LIBSRCS:.f=.o)
	$(AR) $(ARFLAGS) $(M820LIB) $(?:.f=.o)
	$(RANLIB) $(M820LIB)

sources: $(SRCS) $(LIBSRCS)
$(SRCS):
	$(GET) $(REL) $@
$(LIBSRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m820.f
	cat $(SRCS) > ../source/m820.f
	rm -f ../source/src820.f
	cat $(LIBSRCS) > ../source/src820.f

print: link
	lpr ../source/m820.f
	lpr ../source/src820.f

clean:
	rm -f *.o core


