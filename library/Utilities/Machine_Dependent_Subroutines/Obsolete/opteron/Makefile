#   @(#)Makefile	5.2   4/18/95
#
#  Makefile for the machine-dependent portion of MesaLib
#
UTIL_LIB = ../../util_lib_f90.a
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

.f.o:
	$(F90) $(F90FLAGS) $<

all: $(UTIL_LIB)
$(UTIL_LIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(UTIL_LIB) $(?:.f=.o)
	$(RANLIB) $(UTIL_LIB)
sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@

clean:
	rm -f *.o
