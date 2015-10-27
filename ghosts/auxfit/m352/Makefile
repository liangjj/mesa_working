#      %W%   %G%
#
#      Makefile for M352
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
	../xm352/centr3.f \
	../xm352/fixone.f \
	../xm352/fmoft.f \
        ../xm352/form3.f \
	../xm352/gena0int.f\
	../xm352/genacint.f\
	../xm352/loops.f \
	../xm352/m352.f \
	../xm352/overlp.f \
	../xm352/pm352.f \
	../xm352/put3.f \
	../xm352/setcc.f \
	../xm352/sint3.f \
	../xm352/sinvrt.f \
	../xm352/trans3.f \
	../xm352/twoegoof.f\
	../xm352/twoel.f 


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m352

$(BINDIR)/m352: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/m352

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m352.f
	cat $(SRCS) > ../source/m352.f

print: link
	lpr ../source/m352.f

clean:
	rm -f *.o core


