#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M824
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
	../m824/cdens.f \
	../m824/chknrm.f \
	../m824/derfoc.f \
	../m824/dfock.f \
	../m824/fab.f \
	../m824/fmdab.f \
	../m824/lagran.f \
	../m824/lagrng.f \
	../m824/m824.f \
	../m824/mn1010.f \
	../m824/mn333.f \
	../m824/mn814.f \
	../m824/nucrep.f \
	../m824/pm824.f \
	../m824/sder.f \
	../m824/sort32.f \
	../m824/sort64.f \
	../m824/trnd1e.f \
	../m824/trnd2e.f \
	../m824/trnds.f \
	../m824/trnlag.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m824

$(BINDIR)/m824: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m824

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m824.f
	cat $(SRCS) > ../source/m824.f

print: link
	lpr ../source/m824.f

clean:
	rm -f *.o core
