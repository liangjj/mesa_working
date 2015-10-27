#   @(#)Makefile	1.4   8/2/91
#
# Makefile for m7005
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
	../m7005/abcoef.f \
	../m7005/asymn.f \
	../m7005/asymp.f \
	../m7005/c0sqfn.f \
	../m7005/clfun.f \
	../m7005/coulnt.f \
	../m7005/dlfun.f \
	../m7005/drver.f \
	../m7005/extphs.f \
	../m7005/gencou.f \
	../m7005/iregxp.f \
	../m7005/iregxn.f \
	../m7005/locate.f \
	../m7005/lrcoef.f \
	../m7005/m7005.f \
	../m7005/mkx.f \
	../m7005/mkgrd.f \
	../m7005/mkxexp.f \
	../m7005/plfun.f \
	../m7005/potntl.f \
	../m7005/qlplfn.f \
	../m7005/regexn.f \
	../m7005/regexp.f \
	../m7005/rlfun.f \
	../m7005/splinr.f \
	../m7005/splint.f \
	../m7005/spltst.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7005

$(BINDIR)/m7005: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m7005




