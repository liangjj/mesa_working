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
TCGMSGLIB = /usr4/russo/src/pnl_stuff/tcgmsg/libtcgmsg.$(OPSYS).a
MDLIB =
GET = sccs get
REL = 

SRCS = \
	../xm715/d2expn.f \
	../xm715/derexp.f \
	../xm715/dprim.f \
	../xm715/driver.f \
	../xm715/fd2ant.f \
	../xm715/fd2bnt.f \
	../xm715/fmdint.f \
	../xm715/get1dm.f \
	../xm715/getalp.f \
	../xm715/gj2pdm.f \
	../xm715/ldexp.f \
	../xm715/lebedev.f \
	../xm715/m715.f \
	../xm715/makd2g.f \
	../xm715/makd2h.f \
	../xm715/makd2i.f \
	../xm715/makdg.f \
	../xm715/makdh.f \
	../xm715/makdi.f \
	../xm715/md2ghi.f \
	../xm715/mkdghi.f \
	../xm715/movd2e.f \
	../xm715/movder.f \
	../xm715/nxtask.f \
	../xm715/pbasis.f \
	../xm715/plnkerr.f \
	../xm715/pm715.f \
	../xm715/prefac.f \
	../xm715/prims.f \
	../xm715/prntld.f \
	../xm715/redund.f \
	../xm715/sg1.f \
	../xm715/sizer.f \
	../xm715/sphere.f \
	../xm715/vfmint.f \
	../xm715/vmkghi.f 

.f.o:
	$(FC) $(FFLAGS) $<


all: $(BINDIR)/m715 $(BINDIR)/m715.slave

../xm715/getscm.small.o: ../xm715/getscm.small.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) $<

../xm715/getscm.big.o: ../xm715/getscm.big.$(OPSYS).f
	$(FC) -o $@ $(FFLAGS) $<

$(BINDIR)/m715: $(SRCS:.f=.o) $(MESALIB) ../xm715/getscm.big.o
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../xm715/getscm.big.o $(MESALIB) $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m715

$(BINDIR)/m715.slave: $(SRCS:.f=.o) $(MESALIB) ../xm715/getscm.small.o
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ../xm715/getscm.small.o $(MESALIB) $(MDLIB) $(TCGMSGLIB) -o $(BINDIR)/m715.slave

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/xm715.f
	cat $(SRCS) > ../source/xm715.f

print: link
	lpr ../source/xm715.f

clean:
	rm -f *.o core


