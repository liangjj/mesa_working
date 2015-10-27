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
	../ultrix/chain.f \
	../ultrix/comand.f \
	../ultrix/dattim.f \
	../ultrix/getscm.f \
	../ultrix/iadtwp.f \
	../ultrix/intowp.f \
	../ultrix/iocls.f \
	../ultrix/ioget.f \
	../ultrix/iogetc.f \
	../ultrix/ioinq.f \
	../ultrix/ioopn.f \
	../ultrix/ioput.f \
	../ultrix/ioputc.f \
	../ultrix/iorm.f \
	../ultrix/itobyt.f \
	../ultrix/lnkerr.f \
	../ultrix/lxopen.f \
	../ultrix/nchrpw.f \
	../ultrix/shiftl.f \
	../ultrix/shiftr.f \
	../ultrix/stderr.f \
	../ultrix/timing.f \
	../ultrix/tsumry.f \
	../ultrix/unqfil.f \
	../ultrix/usrnam.f \
	../ultrix/versn.f \
	../ultrix/wpadti.f \
	../ultrix/wptbyt.f \
	../ultrix/wptoin.f 

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
	rm -f ../../../source/mdutil-ultrix.f
	cat $(SRCS) > ../../../source/mdutil-ultrix.f

print: link
	lpr ../../../source/mdutil-ultrix.f

clean:
	rm -f *.o
	rm -f core
