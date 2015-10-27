#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6002
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =
BINDIR = $(MESA_BIN)
MESALIB = $(MESA_LIB)/mesalib.a
MDLIB =$(MDLIB)
GET = sccs get
REL = 
GET = sccs get
SRCS = \
	../m6002/c2v.f \
	../m6002/d2h.f \
	../m6002/dmatin.f \
	../m6002/fmch.f \
	../m6002/generf.f \
	../m6002/gfunct.f \
	../m6002/lgndre.f \
	../m6002/nlmxyz.f \
	../m6002/prntv.f \
	../m6002/rdbsis.f \
	../m6002/rdmom.f \
	../m6002/stuff0.f \
	../m6002/stuff1.f \
	../m6002/stuff2.f \
	../m6002/stuff3.f \
	../m6002/stuff4.f \
	../m6002/stuff5.f \
	../m6002/stuff6.f \
	../m6002/vecin.f \
	../m6002/vints.f \
	../m6002/vlr.f \
	../m6002/vstat.f \
	../m6002/wrtbas.f \
	../m6002/wrtpot.f

.f.o:
	$(FC) $(FFLAGS) $<


all: $(BINDIR)/m6002


$(BINDIR)/m395: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m395

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m395.f
	cat $(SRCS) > ../source/m395.f

print: link
	lpr ../source/m395.f

clean:
	rm -f *.o core
