#   @(#)Makefile	1.5   8/3/91
#
# Makefile for M914
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
M914LIB = ../library/m914lib.a
GET = sccs get
REL = 

SRCS = \
	../m914/m914.f \
	../m914/loopdat.f \
	../m914/pm914.f

LIBSRCS = \
 	../m914/fdiag.f \
	../m914/fold.f \
	../m914/foldw.f \
	../m914/getdrt.f \
	../m914/inhmat.f \
	../m914/init4x.f \
	../m914/initdg.f \
	../m914/initex.f \
	../m914/initlp.f \
	../m914/initpg.f \
	../m914/initst.f \
	../m914/initx2.f \
	../m914/loopdat.f \
	../m914/loopyd.f \
	../m914/mn901.f \
	../m914/roots.f \
	../m914/square.f \
	../m914/squarw.f \
	../m914/start.f 

.f.o:
	$(FC) $(FFLAGS)  $<

all: $(BINDIR)/m914
$(BINDIR)/m914: $(SRCS:.f=.o) $(M914LIB) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(M914LIB) $(MESALIB) \
                -o $(BINDIR)/m914
$(M914LIB): $(LIBSRCS:.f=.o)
	$(AR) $(ARFLAGS) $(M914LIB) $(?:.f=.o)
	$(RANLIB) $(M914LIB)

sources: $(SRCS) $(LIBSRCS)
$(SRCS):
	$(GET) $(REL) $@
$(LIBSRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m914.f
	cat $(SRCS) > ../source/m914.f
	rm -f ../source/src914.f
	cat $(LIBSRCS) > ../source/src914.f

print: link
	lpr ../source/m914.f
	lpr ../source/src914.f

clean:
	rm -f *.o core
