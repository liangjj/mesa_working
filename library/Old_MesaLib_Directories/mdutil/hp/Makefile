#   @(#)Makefile	5.1   11/6/94
#
#  Makefile for the machine-dependent portion of MesaLib
#
MESALIB = ../../../library/mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = f77
FFLAGS = -c
GET = sccs get
REL = 

SRCS = \
	../hp/and.f \
	../hp/chain.f \
	../hp/comand.f \
	../hp/dattim.f \
	../hp/getscm.f \
	../hp/iadtwp.f \
	../hp/intowp.f \
	../hp/iocls.f \
	../hp/ioget.f \
	../hp/iogetc.f \
	../hp/ioinq.f \
	../hp/ioopn.f \
	../hp/ioput.f \
	../hp/ioputc.f \
	../hp/iorm.f \
	../hp/itobyt.f \
	../hp/lnkerr.f \
	../hp/lxopen.f \
	../hp/nchrpw.f \
	../hp/or.f \
	../hp/shiftl.f \
	../hp/shiftr.f \
	../hp/sizes.f \
	../hp/stderr.f \
	../hp/timing.f \
	../hp/tsumry.f \
	../hp/unqfil.f \
	../hp/usrnam.f \
	../hp/versn.f \
	../hp/wpadti.f \
	../hp/wptbyt.f \
	../hp/wptoin.f \
	../hp/xor.f

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
	rm -f ../../../source/mdutil-hp.f
	cat $(SRCS) > ../../../source/mdutil-hp.f

print: link
	lpr ../../../source/mdutil-hp.f

clean:
	rm -f *.o
	rm -f core
