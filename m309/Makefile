#      @(#)Makefile	5.1   11/6/94
#
#      Makefile for M309
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
	../m309/fixone.f \
        ../m309/fmonel.f \
        ../m309/fmt.f \
        ../m309/lpints.f \
	../m309/m309.f \
	../m309/nucrep.f \
	../m309/oneint.f \
	../m309/pm309.f \
	../m309/tints.f \
	../m309/vints.f  

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m309

$(BINDIR)/m309: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m309

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m309.f
	cat $(SRCS) > ../source/m309.f

print: link
	lpr ../source/m309.f

clean:
	rm -f *.o core


