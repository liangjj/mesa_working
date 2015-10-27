#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M903
#
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB =
M903LIB = ../library/m903lib.a
GET = sccs get
REL = 

SRCS = \
	../m903/m903.f \
	../m903/pm903.f

LIBSRCS = \
	../m903/bottom.f \
	../m903/ci.f \
	../m903/dm.f \
	../m903/dndiag.f \
	../m903/dnloop.f \
	../m903/dnwalk.f \
	../m903/ecdn.f \
	../m903/ectb.f \
	../m903/ecup.f \
	../m903/fixdm.f \
	../m903/fmdm.f \
	../m903/fmints.f \
	../m903/fmtdm.f \
	../m903/gdrt.f \
	../m903/hc.f \
	../m903/intad.f \
	../m903/mn903.f \
	../m903/multhc.f \
	../m903/onedm.f \
	../m903/oneel.f \
	../m903/onetdm.f \
	../m903/reord.f \
	../m903/rfmint.f \
	../m903/sdn.f \
	../m903/stb.f \
	../m903/sup.f \
	../m903/symdm.f \
	../m903/tdm.f \
	../m903/top.f \
	../m903/twodm.f \
	../m903/twoel.f \
	../m903/twotdm.f \
	../m903/ufmint.f \
	../m903/updrt.f \
	../m903/uploop.f \
	../m903/upwalk.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m903
lib: $(M903LIB)
$(BINDIR)/m903: $(SRCS:.f=.o) $(M903LIB) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(M903LIB) $(MESALIB) \
                $(MDLIB) -o $(BINDIR)/m903
$(M903LIB): $(LIBSRCS:.f=.o)
	$(AR) $(ARFLAGS) $(M903LIB) $(?:.f=.o)
	$(RANLIB) $(M903LIB)

sources: $(SRCS) $(LIBSRCS)
$(SRCS):
	$(GET) $(REL) $@
$(LIBSRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m903.f
	cat $(SRCS) > ../source/m903.f
	rm -f ../source/src903.f
	cat $(LIBSRCS) > ../source/src903.f

print: link
	lpr ../source/m903.f
	lpr ../source/src903.f

clean:
	rm -f *.o core
