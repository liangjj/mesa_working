#   @(#)Makefile	5.2   4/18/95
#
#  Makefile for the machine-dependent portion of MesaLib
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR = $(UTILITY)
include $(TOP_DIR)/Makefile.inc
GET = sccs get
REL = 

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
OBJECTS = \
	$(SRCS:.f=.o) 

all: $(MACHLIB)
$(MACHLIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $(MACHLIB) $(OBJECTS)
	$(RANLIB) $(MACHLIB)
	$(AR) $(ARFLAGS) $(UTILITY_LIB) $(OBJECTS)
	$(RANLIB) $(UTILITY_LIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@
clean:
	rm -f *.o
