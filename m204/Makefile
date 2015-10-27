#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M204
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB =
GET = sccs get
REL = 

SRCS = \
	../m204/eigint.f \
	../m204/freq.f \
	../m204/frmain.f \
	../m204/m204.f \
	../m204/pm204.f \
	../m204/readfx.f \
	../m204/vibfrc.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m204

$(BINDIR)/m204: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m204

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m204.f
	cat $(SRCS) > ../source/m204.f

print: link
	lpr ../source/m204.f

clean:
	rm -f *.o core


