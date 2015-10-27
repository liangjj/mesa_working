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
	../ctss/chain.f \
	../ctss/comand.f \
	../ctss/dattim.f \
	../ctss/getscm.f \
	../ctss/iadtwp.f \
	../ctss/intowp.f \
	../ctss/iocls.f \
	../ctss/ioget.f \
	../ctss/iogetc.f \
	../ctss/ioinq.f \
	../ctss/ioopn.f \
	../ctss/ioput.f \
	../ctss/ioputc.f \
	../ctss/iorm.f \
	../ctss/iounq.f \
	../ctss/itobyt.f \
	../ctss/lnkerr.f \
	../ctss/lxopen.f \
	../ctss/nchrpw.f \
	../ctss/tsumry.f \
	../ctss/unqfil.f \
	../ctss/usrnam.f \
	../ctss/versn.f \
	../ctss/wpadti.f \
	../ctss/wptbyt.f \
	../ctss/wptoin.f 

#all: $(MESALIB)
#$(MESALIB): $(SRCS:.f=.o)
#	$(AR) $(ARFLAGS) $(MESALIB) $(?:.f=.o)
#	$(RANLIB) $(MESALIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@

link: sources
	rm -f ../../../source/mdutil-ctss.f
	cat $(SRCS) > ../../../source/mdutil-ctss.f
	rm -f ../../../source/mesaccl
	cp MESA.CCL ../../../source/mesaccl

print: link
	lpr ../../../source/mdutil-ctss.f

clean:
	rm -f *.o
	rm -f core
