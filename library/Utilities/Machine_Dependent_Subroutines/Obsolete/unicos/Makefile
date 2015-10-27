# @(#)Makefile	5.1 11/6/94 
#
#  Makefile for the machine-dependent portion of MesaLib
#
MESALIB = ../../../library/mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = f77
FFLAGS = -c -c
GET = sccs get
REL = 

SRCS = \
	../unicos/chain.f \
	../unicos/comand.f \
	../unicos/dattim.f \
	../unicos/getscm.f \
	../unicos/iadtwp.f \
	../unicos/intowp.f \
	../unicos/iocls.f \
	../unicos/ioget.f \
	../unicos/iogetc.f \
	../unicos/ioinq.f \
	../unicos/ioopn.f \
	../unicos/ioput.f \
	../unicos/ioputc.f \
	../unicos/iorm.f \
	../unicos/iounq.f \
	../unicos/itobyt.f \
	../unicos/lnkerr.f \
	../unicos/lxopen.f \
	../unicos/nchrpw.f \
	../unicos/stderr.f \
	../unicos/timing.f \
	../unicos/tsumry.f \
	../unicos/unqfil.f \
	../unicos/usrnam.f \
	../unicos/versn.f \
	../unicos/wpadti.f \
	../unicos/wptbyt.f \
	../unicos/wptoin.f 

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
	rm -f ../../../source/mdutil-unicos.f
	cat $(SRCS) > ../../../source/mdutil-unicos.f

print: link
	lpr ../../../source/mdutil-unicos.f

clean:
	rm -f *.o
	rm -f core
