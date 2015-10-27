#   %W%   %G%
#
#  Makefile for the machine-dependent portion of MesaLib
#
MESALIB = ../../../library/mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = f77
FFLAGS = -c -g
GET = sccs get
REL = 

SRCS = \
	../nltss/chain.f \
	../nltss/charges.f \
	../nltss/comand.f \
	../nltss/dattim.f \
	../nltss/eqvfix.f \
	../nltss/getscm.f \
	../nltss/iadtwp.f \
	../nltss/intowp.f \
	../nltss/iocls.f \
	../nltss/ioget.f \
	../nltss/iogetc.f \
	../nltss/ioinq.f \
	../nltss/ioopn.f \
	../nltss/ioput.f \
	../nltss/ioputc.f \
	../nltss/iorm.f \
	../nltss/iounq.f \
	../nltss/itobyt.f \
	../nltss/lnkerr.f \
	../nltss/lxopen.f \
	../nltss/nchrpw.f \
	../nltss/timing.f \
	../nltss/tsumry.f \
	../nltss/unique.f \
	../nltss/unqfil.f \
	../nltss/usrnam.f \
	../nltss/versn.f \
	../nltss/wpadti.f \
	../nltss/wptbyt.f \
	../nltss/wptoin.f 

#all: $(MESALIB)
#$(MESALIB): $(SRCS:.f=.o)
#	$(AR) $(ARFLAGS) $(MESALIB) $(?:.f=.o)
#	$(RANLIB) $(MESALIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@

link: sources
	rm -f ../../../source/mdutil-nltss.f
	cat $(SRCS) > ../../../source/mdutil-nltss.f

print: link
	lpr ../../../source/mdutil-nltss.f

clean:
	rm -f *.o
	rm -f core
