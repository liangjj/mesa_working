#   @(#)Makefile	5.2   4/18/95
#
#  Makefile for the machine-dependent portion of MesaLib
#
MESALIB = ../mesalib.a
AR = -ar
ARFLAGS = rvu
GET = sccs get
REL = 
RANLIB = ranlib

SRCS = \
	asinh.f \
	chain.f \
	comand.f \
	d1mach.f \
	dattim.f \
	dcsevl.f \
	getscm.f \
	getmem.f \
	iadtwp.f \
	initds.f \
	intowp.f \
	iocls.f \
	ioget.f \
	iogetc.f \
	ioinq.f \
	ioopn.f \
	ioput.f \
	ioputc.f \
	iorm.f \
	itobyt.f \
	lnkerr.f \
	lxopen.f \
	nchrpw.f \
	shiftl.f \
	shiftr.f \
	stderr.f \
	timing.f \
	tsumry.f \
	unqfil.f \
	usrnam.f \
	versn.f \
	wind.f \
	wpadti.f \
	wptbyt.f \
	wptoin.f 

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
