#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M830
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
M830LIB = ../library/m830lib.a
GET = sccs get
REL = 

SRCS = \
	../m830/m830.f \
	../m830/pm830.f

LIBSRCS = \
	../m830/canon.f \
	../m830/mn830.f \
	../m830/tocan.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m830
lib: $(M830LIB)
$(BINDIR)/m830: $(SRCS:.f=.o) $(M830LIB) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(M830LIB) $(MESALIB) \
                $(MDLIB) -o $(BINDIR)/m830
$(M830LIB): $(LIBSRCS:.f=.o)
	$(AR) $(ARFLAGS) $(M830LIB) $(?:.f=.o)
	$(RANLIB) $(M830LIB)

sources: $(SRCS) $(LIBSRCS)
$(SRCS):
	$(GET) $(REL) $@
$(LIBSRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m830.f
	cat $(SRCS) > ../source/m830.f
	rm -f ../source/src830.f
	cat $(LIBSRCS) > ../source/src830.f

print: link
	lpr ../source/m830.f
	lpr ../source/src830.f

clean:
	rm -f *.o core


