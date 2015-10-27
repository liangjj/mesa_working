#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M801
#

LD = $(FC)
LDFLAGS = 

BINDIR = $(MESA_BIN)
MESALIB = $(MESA_LIB)/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../m801/bkspac.f \
	../m801/bloksd.f \
	../m801/constd.f \
	../m801/drtcdd.f \
	../m801/drtind.f \
	../m801/getchr.f \
	../m801/getcnt.f \
	../m801/getcod.f \
	../m801/getkey.f \
	../m801/getlin.f \
	../m801/getspe.f \
	../m801/getsym.f \
	../m801/intgrl.f \
	../m801/m801.f \
	../m801/nextc.f \
	../m801/paldus.f \
	../m801/print.f \
	../m801/putdrt.f \
	../m801/ref.f \
	../m801/renumb.f \
	../m801/reordr.f \
	../m801/scan.f \
	../m801/skipbl.f \
	../m801/squeez.f \
	../m801/tapesd.f \
	../m801/weight.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m801

$(BINDIR)/m801: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m801

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m801.f
	cat $(SRCS) > ../source/m801.f

print: link
	lpr ../source/m801.f

clean:
	rm -f *.o core



