#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M902
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
MDLIB = $(MDLIB)
M902LIB = $(MESA_LIB)/m902lib.a
GET = sccs get
REL = 

SRCS = \
	../m902/m902.f \
	../m902/pm902.f

LIBSRCS = \
	../m902/david.f \
	../m902/dm02.f \
	../m902/dmlps.f \
	../m902/eig902.f \
	../m902/exchng.f \
	../m902/fmhpq.f \
	../m902/fmh.f \
	../m902/guess2.f \
	../m902/hamilt.f \
	../m902/hc02.f \
	../m902/hclps.f \
	../m902/hmult.f \
	../m902/incore.f \
	../m902/loops.f \
	../m902/mn902.f \
	../m902/tdm02.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m902
lib: $(M902LIB)
$(BINDIR)/m902: $(SRCS:.f=.o) $(M902LIB) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(M902LIB) $(MESALIB) \
                $(MDLIB) -o $(BINDIR)/m902
$(M902LIB): $(LIBSRCS:.f=.o)
	$(AR) $(ARFLAGS) $(M902LIB) $(?:.f=.o)
	$(RANLIB) $(M902LIB)

sources: $(SRCS) $(LIBSRCS)
$(SRCS):
	$(GET) $(REL) $@
$(LIBSRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m902.f
	cat $(SRCS) > ../source/m902.f
	rm -f ../source/src902.f
	cat $(LIBSRCS) > ../source/src902.f

print: link
	lpr ../source/m902.f
	lpr ../source/src902.f

clean:
	rm -f *.o core
