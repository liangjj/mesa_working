#      @(#)Makefile	5.1   11/6/94
#
#      Makefile for M305
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
	../m305/fixone.f \
        ../m305/lpints.f \
	../m305/ltab.f \
	../m305/m305.f \
	../m305/nucrep.f \
	../m305/oneint.f \
	../m305/pm305.f \
	../m305/tints.f \
	../m305/vints.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m305

$(BINDIR)/m305: $(SRCS:.f=.o) $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m305

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m305.f
	cat $(SRCS) > ../source/m305.f

print: link
	lpr ../source/m305.f

clean:
	rm -f *.o core


