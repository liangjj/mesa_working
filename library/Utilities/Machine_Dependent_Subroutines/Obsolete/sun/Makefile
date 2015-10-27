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
	../sun/asinh.f\
	../sun/chain.f \
	../sun/comand.f \
	../sun/d1mach.f \
	../sun/dattim.f \
	../sun/dcsevl.f \
	../sun/getscm.f \
	../sun/iadtwp.f \
	../sun/initds.f \
	../sun/intowp.f \
	../sun/iocls.f \
	../sun/ioget.f \
	../sun/iogetc.f \
	../sun/ioinq.f \
	../sun/ioopn.f \
	../sun/ioput.f \
	../sun/ioputc.f \
	../sun/iorm.f \
	../sun/itobyt.f \
	../sun/lnkerr.f \
	../sun/lxopen.f \
	../sun/nchrpw.f \
	../sun/shiftl.f \
	../sun/shiftr.f \
	../sun/stderr.f \
	../sun/timing.f \
	../sun/tsumry.f \
	../sun/unqfil.f \
	../sun/usrnam.f \
	../sun/versn.f \
	../sun/wpadti.f \
	../sun/wptbyt.f \
	../sun/wptoin.f 

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
	rm -f ../../../source/mdutil-sun.f
	cat $(SRCS) > ../../../source/mdutil-sun.f

print: link
	lpr ../../../source/mdutil-sun.f

clean:
	rm -f *.o
	rm -f core
