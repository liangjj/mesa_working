#   @(#)Makefile	1.4   8/1/91
#
# Makefile for M910
#
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
FC = f77
FFLAGS = -c -O4

LD = f77
LDFLAGS =

BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
LIBDXML = ../library/libdxml.a
DVDLIB = ../m910/dvdson/dvdlib.a
M910LIB = ../library/m910lib.a
GET = sccs get
REL = 

SRCS = \
	../m910/m910.f

LIBSRCS = \
	../m910/david.f \
	../m910/diagh.f \
	../m910/eig910.f \
	../m910/exchng.f \
	../m910/guess2.f \
	../m910/hamilt.f \
	../m910/hmult.f \
	../m910/incore.f \
	../m910/loops.f \
	../m910/mn910.f \
	../m910/rdham.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m910
lib: $(M910LIB)
$(BINDIR)/m910: $(SRCS:.f=.o) $(M910LIB) $(DVDLIB) $(MESCTLIB) $(LIBDXML)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(M910LIB) $(DVDLIB) $(MESCTLIB) \
	                               $(LIBDXML) -o $(BINDIR)/m910
$(M910LIB): $(LIBSRCS:.f=.o)
	$(AR) $(ARFLAGS) $(M910LIB) $(?:.f=.o)
	$(RANLIB) $(M910LIB)

sources: $(SRCS) $(LIBSRCS)
$(SRCS):
	$(GET) $(REL) $@
$(LIBSRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m910.f
	cat $(SRCS) > ../source/m910.f
	rm -f ../source/src910.f
	cat $(LIBSRCS) > ../source/src910.f

print: link
	lpr ../source/m910.f
	lpr ../source/src910.f

clean:
	rm -f *.o core
















