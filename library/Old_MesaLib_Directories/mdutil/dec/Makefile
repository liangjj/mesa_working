#   @(#)Makefile	5.2   4/18/95
#
#  Makefile for the machine-dependent portion of MesaLib
#
MESALIB = ../../../library/mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = $(FC)
FFLAGS = $(FFLAGS)
GET = sccs get
REL = 

SRCS = \
	../dec/asinh.f \
	../dec/chain.f \
	../dec/comand.f \
	../dec/d1mach.f \
	../dec/dattim.f \
	../dec/dcsevl.f \
	../dec/getscm.f \
	../dec/getmem.f \
	../dec/iadtwp.f \
	../dec/initds.f \
	../dec/intowp.f \
	../dec/iocls.f \
	../dec/ioget.f \
	../dec/iogetc.f \
	../dec/ioinq.f \
	../dec/ioopn.f \
	../dec/ioput.f \
	../dec/ioputc.f \
	../dec/iorm.f \
	../dec/itobyt.f \
	../dec/lnkerr.f \
	../dec/lxopen.f \
	../dec/nchrpw.f \
	../dec/shiftl.f \
	../dec/shiftr.f \
	../dec/stderr.f \
	../dec/timing.f \
	../dec/tsumry.f \
	../dec/unqfil.f \
	../dec/usrnam.f \
	../dec/versn.f \
	../dec/wind.f \
	../dec/wpadti.f \
	../dec/wptbyt.f \
	../dec/wptoin.f 

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
	rm -f ../../../source/mdutil-dec.f
	cat $(SRCS) > ../../../source/mdutil-dec.f

print: link
	lpr ../../../source/mdutil-dec.f

clean:
	rm -f *.o
	rm -f core
