#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M250
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
	../m250/eigint.f \
	../m250/freq.f \
	../m250/frmain.f \
	../m250/getgrd.f \
	../m250/m250.f \
	../m250/pm250.f \
	../m250/project.f \
	../m250/rphdrv.f \
	../m250/stepp.f \
	../m250/vibfrc.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m250

$(BINDIR)/m250: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m250

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m250.f
	cat $(SRCS) > ../source/m250.f

print: link
	lpr ../source/m250.f

clean:
	rm -f *.o core



