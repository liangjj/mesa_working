#   @(#)Makefile	1.3   8/1/91
#
# Makefile for Resfit
#
FC = fort
FFLAGS = -c -O4

LD = fort
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB = /usr/lib/libcxml.a
GET = sccs get
REL = 

SRCS = \
	../m6025/bnds.f \
	../m6025/condens.f \
	../m6025/datain.f \
	../m6025/der.f \
	../m6025/expad.f \
	../m6025/fcn.f \
	../m6025/gatan.f \
	../m6025/glss.f  \
	../m6025/gtan.f \
	../m6025/lstsq.f \
	../m6025/mmax.f \
	../m6025/normal.f \
	../m6025/opac.f \
	../m6025/phasit.f \
	../m6025/phig.f \
	../m6025/pnear.f \
	../m6025/report.f \
	../m6025/resfit.f \
	../m6025/rownorm.f \
	../m6025/srch.f \
	../m6025/stest.f \
	../m6025/swriti.f \
	../m6025/swritd.f \
	../m6025/timefs.f \
	../m6025/unpak.f \
	../m6025/vas.f


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6025

$(BINDIR)/m6025: $(SRCS:.f=.o) $(MESALIB) $(MDLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m6025

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m6025.f
	cat $(SRCS) > ../source/m6025.f

print: link
	lpr ../source/m101.f

clean:
	rm -f *.o core
