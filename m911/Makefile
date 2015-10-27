#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M911
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
	../m911/dot.f \
	../m911/fdiag.f \
	../m911/getdrt.f \
	../m911/inimat.f \
	../m911/init4x.f \
	../m911/initdg.f \
	../m911/initex.f \
	../m911/initlp.f \
	../m911/initnx.f \
	../m911/initpg.f \
	../m911/initst.f \
	../m911/initx2.f \
	../m911/loopyd.f \
	../m911/lupdat.f \
	../m911/m911.f \
	../m911/pm911.f \
	../m911/putdm.f \
	../m911/square.f \
	../m911/squarw.f \
	../m911/start.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m911

$(BINDIR)/m911: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m911

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m911.f
	cat $(SRCS) > ../source/m911.f

print: link
	lpr ../source/m911.f

clean:
	rm -f *.o core



