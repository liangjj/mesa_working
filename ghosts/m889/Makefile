#   @(#)Makefile	1.1   11/30/90
#
# Makefile for M889
#
FC = f77
FFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../m889/chknrm.f \
	../m889/fock.f \
	../m889/getdrt.f \
	../m889/m889.f \
	../m889/pm811.f \
	../m889/reordr.f \
	../m889/trn1e.f \
	../m889/trn2e.f \
	../m889/vecout.f 

.f.o:
	$(FC) $(FFLAGS) -c $<

all: $(BINDIR)/m889

$(BINDIR)/m889: $(SRCS:.f=.o) $(MESALIB)
	$(FC) $(FFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/m889

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m889.f
	cat $(SRCS) > ../source/m889.f

print: link
	lpr ../source/m889.f

clean:
	rm -f *.o core



