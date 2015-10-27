#   @(#)Makefile	4.2   7/16/93
#
#  Makefile for the I/O portion of MesaLib
#
MESALIB = ../mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB =  ranlib

FC = f77
FFLAGS = -c
F77FLAGS = $(FFLAGS) +E1 +ppu
GET = sccs get
REL = 

SRCS = \
	../io/ctoi.f \
	../io/dcaptl.f \
	../io/dirsrt.f \
	../io/ffnext.f \
	../io/gettok.f \
	../io/ioabor.f \
	../io/ioclos.f \
	../io/iocls.f \
	../io/iocopy.f \
	../io/iodest.f \
	../io/iodump.f \
	../io/ioeof.f \
	../io/ioexst.f \
	../io/iofile.f \
	../io/ioflno.f \
	../io/ioget.f \
	../io/iogetc.f \
	../io/ioinq.f \
	../io/iolen.f \
	../io/ioopen.f \
	../io/ioopn.f \
	../io/ioput.f \
	../io/ioputc.f \
	../io/ioread.f \
	../io/iorew.f \
	../io/iorm.f \
	../io/iosize.f \
	../io/iosys.f \
	../io/iounit.f \
	../io/iowait.f \
	../io/iowrit.f \
	../io/iowt.f \
	../io/iowtab.f \
	../io/itobyt.f  \
	../io/itoc.f  \
	../io/lnkerr.f  \
	../io/unitwt.f  \
	../io/wptbyt.f	\
	../io/stderr.f	\
	../io/wind.f

.f.o:
	$(FC) -c $(F77FLAGS) $<

all: $(MESALIB)
$(MESALIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(MESALIB) $(?:.f=.o)
	$(RANLIB) $(MESALIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@

link: sources
	rm -f ../../source/io.f
	cat $(SRCS) > ../../source/io.f

print: link
	lpr ../../source/io.f

clean:
	rm -f *.o
	rm -f core

