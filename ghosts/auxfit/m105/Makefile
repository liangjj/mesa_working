#   %W%   %G%
#
# Makefile for M105
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../xm105/baspr1.f \
	../xm105/basprt.f \
	../xm105/ecpprt.f \
	../xm105/genaux.f \
	../xm105/genbas.f \
	../xm105/getbas.f \
	../xm105/grppt.f \
	../xm105/grppt1.f \
	../xm105/m105.f \
	../xm105/maket.f \
	../xm105/maket1.f \
	../xm105/moment.f \
	../xm105/pm105.f \
	../xm105/symptr.f


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/xm105

$(BINDIR)/xm105: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/m105

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/xm105.f
	cat $(SRCS) > ../source/xm105.f

print: link
	lpr ../source/xm105.f

clean:
	rm -f *.o core


