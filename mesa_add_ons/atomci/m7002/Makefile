#  @(#)Makefile	1.8 8/1/91
#
# Makefile for M7002
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =

BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a

GET = sccs get
REL = 

SRCS = \
	../m7002/clebgd.f \
	../m7002/clebs.f \
	../m7002/racah.f \
	../m7002/rdcleb.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7002

$(BINDIR)/m7002: $(SRCS:.f=.o) $(MESCTLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) -o $(BINDIR)/m7002

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m7002.f
	cat $(SRCS) > ../source/m7002.f

print: link
	lpr ../source/m7002.f
clean:
	rm -f *.o core




