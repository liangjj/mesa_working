#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M5000
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
	cham.f \
	chk.f \
	compar.f \
	convt.f \
	diagx.f \
	drslv.f \
	dropfn.f \
	drvlin.f \
	honv.f \
	ihamxt.f \
	lin.f \
	linslv.f \
	mkrhs.f \
	prepfn.f \
	rdham.f \
	rdiag.f \
	setind.f \
	tcart.f \
	vmat.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m5000

$(BINDIR)/m5000: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m5000

