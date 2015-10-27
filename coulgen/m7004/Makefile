#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m7004
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
	../m7004/abcoef.f \
	../m7004/asymn.f \
	../m7004/asymp.f \
	../m7004/c0sqfn.f \
	../m7004/clfun.f \
	../m7004/coulnt.f \
	../m7004/dlfun.f \
	../m7004/drver.f \
	../m7004/extcou.f \
	../m7004/extphs.f \
	../m7004/fndswt.f \
	../m7004/fninpt.f \
	../m7004/gencou.f \
	../m7004/inhomo.f \
	../m7004/iregxp.f \
	../m7004/iregxn.f \
	../m7004/lrcoef.f \
	../m7004/m7004.f \
	../m7004/mkgrd.f \
	../m7004/mkexp.f \
	../m7004/plfun.f \
	../m7004/potntl.f \
	../m7004/qlplfn.f \
	../m7004/regexn.f \
	../m7004/regexp.f \
	../m7004/rlfun.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7004

$(BINDIR)/m7004: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m7004




