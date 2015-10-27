#   @(#)Makefile	5.2   4/18/95
#
#  Makefile for the machine-dependent portion of MesaLib
#
MESALIB = ../../../library/mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c
GET = sccs get
REL = 

SRCS = \
	../sgi/asinh.f \
	../sgi/chain.f \
	../sgi/comand.f \
	../sgi/d1mach.f \
	../sgi/dattim.f \
	../sgi/dcsevl.f \
	../sgi/getscm.f \
	../sgi/iadtwp.f \
	../sgi/initds.f \
	../sgi/intowp.f \
	../sgi/iocls.f \
	../sgi/ioget.f \
	../sgi/iogetc.f \
	../sgi/ioinq.f \
	../sgi/ioopn.f \
	../sgi/ioput.f \
	../sgi/ioputc.f \
	../sgi/iorm.f \
	../sgi/itobyt.f \
	../sgi/lnkerr.f \
	../sgi/lxopen.f \
	../sgi/memory.f \
	../sgi/nchrpw.f \
	../sgi/shiftl.f \
	../sgi/shiftr.f \
	../sgi/stderr.f \
	../sgi/timing.f \
	../sgi/tsumry.f \
	../sgi/unqfil.f \
	../sgi/usrnam.f \
	../sgi/versn.f \
	../sgi/wind.f \
	../sgi/wpadti.f \
	../sgi/wptbyt.f \
	../sgi/wptoin.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(MESALIB)
$(MESALIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(MESALIB) $(?:.f=.o)
	$(RANLIB) $(MESALIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@

link: sources
	rm -f ../../../source/mdutil-sgi.f
	cat $(SRCS) > ../../../source/mdutil-sgi.f

print: link
	lpr ../../../source/mdutil-sgi.f

clean:
	rm -f *.o
	rm -f core
