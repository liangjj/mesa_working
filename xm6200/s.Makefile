h16146
s 00072/00000/00000
d D 1.1 94/02/16 20:34:55 mesa 1 0
c date and time created 94/02/16 20:34:55 by mesa
e
u
U
f e 0
t
T
I 1
# %W% %G%
#
# Makefile for M6200
#
FC = f77
FFLAGS = -c 

LD = f77
LDFLAGS =
BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB = 

GET = sccs get

SRCS = \
	../m6200/aij.f \
	../m6200/ang.f \
	../m6200/cardin.f \
	../m6200/class.f \
	../m6200/gaussq.f \
	../m6200/gamfun.f \
	../m6200/gbslve.f \
	../m6200/gbtql2.f \
	../m6200/getbnd.f \
	../m6200/getcnd.f \
	../m6200/getdnd.f \
	../m6200/lebdev.f \
	../m6200/m6200.f \
	../m6200/mkgr.f \
	../m6200/mkvwt.f \
	../m6200/mkwt.f \
	../m6200/mkyunt.f \
	../m6200/necote.f \
	../m6200/pm6200.f \
	../m6200/radquad.f \
	../m6200/rdang.f \
	../m6200/rdlebdv.f \
	../m6200/rdpt.f \
	../m6200/satshl.f \
	../m6200/scalwt.f \
	../m6200/shells.f \
	../m6200/sumncw.f \
	../m6200/voronoi.f \
	../m6200/xyz.f \
	../m6200/xyzw.f \
	../m6200/yukawa.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6200

$(BINDIR)/m6200: $(SRCS:.f=.o) $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m6200

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m6200.f
	cat $(SRCS) > ../source/m6200.f

print: link
	lpr ../source/m6200.f

clean:
	rm -f *.o core



E 1
