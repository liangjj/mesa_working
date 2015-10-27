#    @(#)Makefile	5.1 11/6/94
#
# Makefile for M821
#
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
FC = $(FC)
FFLAGS = $(FFLAGS)

LD = $(FC)
LDFLAGS =

BINDIR = $(MESA_BIN)
MESALIB = $(MESA_LIB)/mesalib.a
MDLIB =
M821LIB = $(MESA_LIB)/m821lib.a
GET = sccs get
REL = 

SRCS = \
	../m821/m821.f \
	../m821/pm821.f

LIBSRCS = \
	../m821/drt821.f \
	../m821/mn821.f \
	../m821/xxint.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m821
lib: $(M821LIB)
$(BINDIR)/m821: $(SRCS:.f=.o) $(M821LIB) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(M821LIB) $(MESALIB) \
                $(MDLIB) -o $(BINDIR)/m821
$(M821LIB): $(LIBSRCS:.f=.o)
	$(AR) $(ARFLAGS) $(M821LIB) $(?:.f=.o)
	$(RANLIB) $(M821LIB)

sources: $(SRCS) $(LIBSRCS)
$(SRCS):
	$(GET) $(REL) $@
$(LIBSRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m821.f
	cat $(SRCS) > ../source/m821.f
	rm -f ../source/src821.f
	cat $(LIBSRCS) > ../source/src821.f

print: link
	lpr ../source/m821.f
	lpr ../source/src821.f

clean:
	rm -f *.o core


