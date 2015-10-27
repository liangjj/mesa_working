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
	../cos/chain.f \
	../cos/comand.f \
	../cos/dattim.f \
	../cos/getscm.f \
	../cos/iadtwp.f \
	../cos/intowp.f \
	../cos/iocls.f \
	../cos/ioget.f \
	../cos/iogetc.f \
	../cos/ioinq.f \
	../cos/ioopn.f \
	../cos/ioput.f \
	../cos/ioputc.f \
	../cos/iorm.f \
	../cos/iounq.f \
	../cos/itobyt.f \
	../cos/lnkerr.f \
	../cos/lxopen.f \
	../cos/nchrpw.f \
	../cos/tsumry.f \
	../cos/unqfil.f \
	../cos/usrnam.f \
	../cos/versn.f \
	../cos/wpadti.f \
	../cos/wptbyt.f \
	../cos/wptoin.f 

#all: $(MESALIB)
#$(MESALIB): $(SRCS:.f=.o)
#	$(AR) $(ARFLAGS) $(MESALIB) $(?:.f=.o)
#	$(RANLIB) $(MESALIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@

link: sources
	rm -f ../../../source/mdutil-cos.f
	cat $(SRCS) > ../../../source/mdutil-cos.f

print: link
	lpr ../../../source/mdutil-cos.f

clean:
	rm -f *.o
	rm -f core
