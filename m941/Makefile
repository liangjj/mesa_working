#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M901
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
M901LIB = ../library/m901lib.a
GET = sccs get
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

sources: $(SRCS) $(LIBSRCS)
$(SRCS):
	$(GET) $(REL) $@
$(LIBSRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m901.f
	cat $(SRCS) > ../source/m901.f
	rm -f ../source/src901.f
	cat $(LIBSRCS) > ../source/src901.f

print: link
	lpr ../source/m901.f
	lpr ../source/src901.f

clean:
	rm -f *.o core
