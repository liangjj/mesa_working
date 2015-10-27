#   @(#)Makefile	1.4   8/1/91
#
# Makefile for M916
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
M916LIB = ../library/m916lib.a
GET = sccs get
REL = 

SRCS = \
	../m916/m916.f \
	../m916/loopdat.f \
	../m916/pm916.f

LIBSRCS = \
 	../m916/fdiag.f \
	../m916/fold.f \
	../m916/foldw.f \
	../m916/getdrt.f \
	../m916/inhmat.f \
	../m916/init4x.f \
	../m916/initdg.f \
	../m916/initex.f \
	../m916/initlp.f \
	../m916/initpg.f \
	../m916/initst.f \
	../m916/initx2.f \
	../m916/loopdat.f \
	../m916/loopyd.f \
	../m916/mn901.f \
	../m916/roots.f \
	../m916/square.f \
	../m916/squarw.f \
	../m916/start.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m916
$(BINDIR)/m916: $(SRCS:.f=.o) $(M916LIB) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(M916LIB) $(MESALIB) \
                -o $(BINDIR)/m916
$(M916LIB): $(LIBSRCS:.f=.o)
	$(AR) $(ARFLAGS) $(M916LIB) $(?:.f=.o)
	$(RANLIB) $(M916LIB)

sources: $(SRCS) $(LIBSRCS)
$(SRCS):
	$(GET) $(REL) $@
$(LIBSRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m916.f
	cat $(SRCS) > ../source/m916.f
	rm -f ../source/src916.f
	cat $(LIBSRCS) > ../source/src916.f

print: link
	lpr ../source/m916.f
	lpr ../source/src916.f

clean:
	rm -f *.o core
