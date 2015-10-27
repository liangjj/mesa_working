#   @(#)Makefile	5.2   11/28/95
#
#  Makefile for the machine-dependent portion of MesaLib
#
MESALIB = ../../../library/mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = f77
FFLAGS = -c
CC = cc
AS = as
CFLAGS = -c 
AFLAGS = 
GET = sccs get
REL = 

SRCS = \
	../aix/chain.f \
	../aix/comand.f \
	../aix/dattim.f \
	../aix/getscm.f \
	../aix/iadtwp.f \
	../aix/intowp.f \
	../aix/iocls.f \
	../aix/ioget.f \
	../aix/iogetc.f \
	../aix/ioinq.f \
	../aix/ioopn.f \
	../aix/ioput.f \
	../aix/ioputc.f \
	../aix/iorm.f \
	../aix/itobyt.f \
	../aix/lnkerr.f \
	../aix/lxopen.f \
	../aix/movdt.f \
	../aix/nchrpw.f \
	../aix/shiftl.f \
	../aix/shiftr.f \
	../aix/stderr.f \
	../aix/timing.f \
	../aix/tsumry.f \
	../aix/unqfil.f \
	../aix/usrnam.f \
	../aix/versn.f \
	../aix/wpadti.f \
	../aix/wptbyt.f \
	../aix/wptoin.f 


CSRCS = \
	../aix/second.c \
	../aix/loc.c

#OBSOLETE AND NOT PORTABLE!  (I.E. not portable to new IBMS!)
ASRCS = \
	../aix/readrtc.s

.f.o:
	$(FC) $(FFLAGS) $<
.c.o:
	$(CC) $(CFLAGS) $<
.s.o:
	$(AS) $(AFLAGS) -o $(ASRCS:.s=.o) $<

all: $(MESALIB)
$(MESALIB): $(SRCS:.f=.o) $(CSRCS:.c=.o) 
	$(AR) $(ARFLAGS) $(MESALIB) $(?:.f=.o)
	$(RANLIB) $(MESALIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@

link: sources
	rm -f ../../../source/mdutil-aix.f
	rm -f ../../../source/mdutil-aix.c
	rm -f ../../../source/mdutil-aix.s
	cat $(SRCS) > ../../../source/mdutil-aix.f
	cat $(CSRCS) > ../../../source/mdutil-aix.c
	cat $(ASRCS) > ../../../source/mdutil-aix.s

print: link
	lpr ../../../source/mdutil-aix.f
	lpr ../../../source/mdutil-aix.c
	lpr ../../../source/mdutil-aix.s

clean:
	rm -f *.o
	rm -f core
