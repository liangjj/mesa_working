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
	../stardent/and.f \
        ../stardent/asinh.f \
	../stardent/chain.f \
	../stardent/comand.f \
        ../stardent/d1mach.f \
	../stardent/dattim.f \
        ../stardent/dcsevl.f \
	../stardent/getscm.f \
	../stardent/iadtwp.f \
        ../stardent/initds.f \
	../stardent/intowp.f \
	../stardent/iocls.f \
	../stardent/ioget.f \
	../stardent/iogetc.f \
	../stardent/ioinq.f \
	../stardent/ioopn.f \
	../stardent/ioput.f \
	../stardent/ioputc.f \
	../stardent/iorm.f \
	../stardent/itobyt.f \
	../stardent/lnkerr.f \
	../stardent/lxopen.f \
	../stardent/nchrpw.f \
	../stardent/or.f \
	../stardent/shiftl.f \
	../stardent/shiftr.f \
	../stardent/stderr.f \
	../stardent/timing.f \
	../stardent/tsumry.f \
	../stardent/unqfil.f \
	../stardent/usrnam.f \
	../stardent/versn.f \
	../stardent/wpadti.f \
	../stardent/wptbyt.f \
	../stardent/wptoin.f  \
	../stardent/xor.f

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
	rm -f ../../../source/mdutil-stardent.f
	cat $(SRCS) > ../../../source/mdutil-stardent.f

print: link
	lpr ../../../source/mdutil-stardent.f

clean:
	rm -f *.o
	rm -f core
