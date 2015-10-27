#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M5010
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
	diagx.f \
	drvlst.f \
	ham.f \
	honv.f \
	lstsq.f \
	mkrhs.f \
	precon.f \
	seidel.f \
	tstovl.f \
        tstsol.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m5010

$(BINDIR)/m5010: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m5010




