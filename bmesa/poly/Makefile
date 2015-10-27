#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6245
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
	davdag.f \
	david.f \
	diag.f \
	diagx.f \
	diaslv.f \
	diff.f \
	exact.f \
	fcoef.f \
	guess.f \
	honv.f \
	lancz.f \
	matx.f \
	mkwts.f \
        mktrid.f \
	modcof.f \
	oldin.f \
	phase.f \
	poly.f \
	seidel.f \
	toorth.f \
	tridag.f \
	trnham.f \
	trnply.f \
	wtfun.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6245

$(BINDIR)/m6245: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6245
