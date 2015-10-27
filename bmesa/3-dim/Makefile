#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6275
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = /newdsk/mesa/bin
MESCTLIB = /newdsk/mesa/library/mesctlib.a
VECLIB = /usr/lib/libvec.a

GET = sccs get

SRCS = \
	chk.f \
	convt.f \
	david.f \
	davdag.f \
	diagx.f \
	diaslv.f \
	filprj.f \
	guess.f \
	hamil.f \
	honv.f \
	rdham.f \
	rmat.f \
	set2x.f \
        set3x.f \
	setind.f \
	surf.f \
	threed.f \
	tmat.f \
	zguess.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6275

$(BINDIR)/m6275: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6275
