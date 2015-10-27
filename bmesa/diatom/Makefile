#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6255
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
	ddiag.f \
	diagx.f \
	diatom.f \
        filham.f \
        honc.f \
	itdiag.f \
        kethe.f \
        ker.f \
        potntl.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6255

$(BINDIR)/m6255: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6255
