#   @(#)Makefile	1.7   8/7/91
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
	../dec/and.f \
	../dec/chain.f \
	../dec/comand.f \
	../dec/dattim.f \
	../dec/getscm.f \
	../dec/iadtwp.f \
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
	../dec/manmem.f \
	../dec/memory.f \
	../dec/nchrpw.f \
	../dec/or.f \
	../dec/shiftl.f \
	../dec/shiftr.f \
	../dec/timing.f \
	../dec/tsumry.f \
	../dec/unqfil.f \
	../dec/usrnam.f \
	../dec/versn.f \
	../dec/wpadti.f \
	../dec/wptbyt.f \
	../dec/wptoin.f \
	../dec/xor.f

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
