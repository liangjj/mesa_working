#   @(#)Makefile	1.1   11/30/90
#
# Makefile for M841
#
FC = f77
FFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../m841/chk841.f \
	../m841/m841.f \
	../m841/mn841.f \
	../m841/pm841.f \
	../m841/toder.f 

.f.o:
	$(FC) $(FFLAGS) -c $<

all: $(BINDIR)/m841

$(BINDIR)/m841: $(SRCS:.f=.o) $(MESALIB)
	$(FC) $(FFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/m841

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m841.f
	cat $(SRCS) > ../source/m841.f

print: link
	lpr ../source/m841.f

clean:
	rm -f *.o core




