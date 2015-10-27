#   @(#)Makefile	1.1   11/30/90
#
# Makefile for M890
#
AR = -ar
RANLIB = ranlib
FC = f77
FFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
M890LIB = ../library/m890lib.a
GET = sccs get
REL = 

SRCS = \
	../m890/m890.f \
	../m890/pm890.f

LIBSRCS = \
	../m890/guglab.f \
	../m890/mn890.f \
	../m890/sriiii.f \
	../m890/sriijj.f \
	../m890/srijij.f \
	../m890/srijkl.f \
	../m890/srone.f \
	../m890/toguga.f 

.f.o:
	$(FC) $(FFLAGS) -c $<

all: $(BINDIR)/m890
$(BINDIR)/m890: $(SRCS:.f=.o) $(M890LIB) $(MESALIB)
	$(FC) $(FFLAGS) $(SRCS:.f=.o) $(M890LIB) $(MESALIB)\
                 -o $(BINDIR)/m890
$(M890LIB): $(LIBSRCS:.f=.o)
	$(AR) rvu $(M890LIB) $(?:.f=.o)
	$(RANLIB) $(M890LIB)

sources: $(SRCS) $(LIBSRCS)
$(SRCS):
	$(GET) $(REL) $@
$(LIBSRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m890.f
	cat $(SRCS) > ../source/m890.f
	rm -f ../source/src890.f
	cat $(LIBSRCS) > ../source/src890.f

print: link
	lpr ../source/m890.f
	lpr ../source/src890.f

clean:
	rm -f *.o core


