#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6250
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
        conrmt.f \
	diag.f \
	egnpse.f \
	fildag.f \
	filham.f \
	filodg.f \
	filpot.f \
	filv.f \
	match.f \
        mkcfun.f \
        potntl.f \
        prep.f \
        renrm.f \
	rmtrx.f \
        smat.f \
        surpsi.f \
	tmat.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6250

$(BINDIR)/m6250: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6250
