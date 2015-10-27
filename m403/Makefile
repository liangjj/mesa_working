#  @(#)Makefile	5.1  11/6/94 
#
# Makefile for M403
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
	../m403/m403.f \
	../m403/ovrlap.f \
	../m403/pm403.f \
	../m403/prjovr.f \
	../m403/putrec.f \
	../m403/rdbas.f \
	../m403/rdmo.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m403

$(BINDIR)/m403: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m403

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m403.f
	cat $(SRCS) > ../source/m403.f

print: link
	lpr ../source/m403.f

clean:
	rm -f *.o core


