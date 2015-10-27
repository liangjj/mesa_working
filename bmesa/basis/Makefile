#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m6237
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../../bin
MESCTLIB = ../../library/mesctlib.a
VECLIB = /usr/lib/libvec.a
GET = sccs get

SRCS = \
	../tridiag/addend.f \
	../tridiag/bndcnd.f \
	../tridiag/better.f \
	../tridiag/condit.f \
	../tridiag/dtafd.f \
	../tridiag/fdiff.f \
	../tridiag/fzero.f \
	../tridiag/gfdiff.f \
        ../tridiag/hdiag.f \
	../tridiag/mkgrd.f \
	../tridiag/nwtprn.f \
	../tridiag/orth.f \
        ../tridiag/potntl.f \
	../tridiag/richdn.f \
	../tridiag/rmteig.f \
	../tridiag/sumncw.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6237

$(BINDIR)/m1205: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6237




