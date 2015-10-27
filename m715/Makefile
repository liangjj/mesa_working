#   %W%   %G%
#
# Makefile for M715
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
TCGMSGLIB = $(MESA_LIBS)/tcgmsg/libtcgmsg.$(OPSYS).a
MDLIB =
GET = sccs get
REL = 

SRCS = \
	../m715/d2expn.f \
	../m715/derexp.f \
	../m715/dprim.f \
	../m715/driver.f \
	../m715/fd2ant.f \
	../m715/fd2bnt.f \
	../m715/fmdint.f \
	../m715/get1dm.f \
	../m715/getalp.f \
	../m715/gj2pdm.f \
	../m715/ldexp.f \
	../m715/lebedev.f \
	../m715/m715.f \
	../m715/makd2g.f \
	../m715/makd2h.f \
	../m715/makd2i.f \
	../m715/makdg.f \
	../m715/makdh.f \
	../m715/makdi.f \
	../m715/md2ghi.f \
	../m715/mkdghi.f \
	../m715/movd2e.f \
	../m715/movder.f \
	../m715/nxtask.f \
	../m715/pbasis.f \
	../m715/plnkerr.f \
	../m715/pm715.f \
	../m715/prefac.f \
	../m715/prims.f \
	../m715/prntld.f \
	../m715/redund.f \
	../m715/sg1.f \
	../m715/sizer.f \
	../m715/sphere.f \
	../m715/vfmint.f \
	../m715/vmkghi.f 

.f.o:
	$(FC) $(FFLAGS) $<


all: $(BINDIR)/m715 $(BINDIR)/m715.slave

../m715/getscm.small.o: ../m715/getscm.small.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) ../m715/getscm.small.$(OPSYS).f

../m715/getscm.big.o: ../m715/getscm.big.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) ../m715/getscm.big.$(OPSYS).f 

$(BINDIR)/m715: $(SRCS:.f=.o) $(MESALIB) ../m715/getscm.big.o
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../m715/getscm.big.o $(MESALIB) $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m715

$(BINDIR)/m715.slave: $(SRCS:.f=.o) $(MESALIB) ../m715/getscm.small.o
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../m715/getscm.small.o $(MESALIB) $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m715.slave

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m715.f
	cat $(SRCS) > ../source/m715.f

print: link
	lpr ../source/m715.f

clean:
	rm -f *.o core


