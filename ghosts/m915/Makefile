#   @(#)Makefile	2.1   10/10/91
#
# Makefile for M915
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
M915LIB = ../library/m915lib.a
GET = sccs get
REL = 

SRCS = \
	../m915/m915.f \
	../m915/loopdat.f \
	../m915/pm915.f

LIBSRCS = \
 	../m915/fdiag.f \
	../m915/fold.f \
	../m915/foldw.f \
	../m915/getdrt.f \
	../m915/inhmat.f \
	../m915/init4x.f \
	../m915/initdg.f \
	../m915/initex.f \
	../m915/initlp.f \
	../m915/initpg.f \
	../m915/initst.f \
	../m915/initx2.f \
	../m915/loopdat.f \
	../m915/loopyd.f \
	../m915/mn901.f \
	../m915/roots.f \
	../m915/square.f \
	../m915/squarw.f \
	../m915/start.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m915
$(BINDIR)/m915: $(SRCS:.f=.o) $(M915LIB) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(M915LIB) $(MESALIB) \
                -o $(BINDIR)/m915
$(M915LIB): $(LIBSRCS:.f=.o)
	$(AR) $(ARFLAGS) $(M915LIB) $(?:.f=.o)
	$(RANLIB) $(M915LIB)

sources: $(SRCS) $(LIBSRCS)
$(SRCS):
	$(GET) $(REL) $@
$(LIBSRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m915.f
	cat $(SRCS) > ../source/m915.f
	rm -f ../source/src915.f
	cat $(LIBSRCS) > ../source/src915.f

print: link
	lpr ../source/m915.f
	lpr ../source/src915.f

clean:
	rm -f *.o core
