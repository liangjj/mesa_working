#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M6050
#
FC = f77
FFLAGS = -c +OP4

LD = f77
LDFLAGS =
BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a
GET = sccs get
SRCS = \
	../m6050/elm.f \
	../m6050/grdprn.f \
	../m6050/ilm.f \
	../m6050/mkwlm.f \
	../m6050/mkthet.f \
	../m6050/nsvlam.f \
	../m6050/opt2e.f \
	../m6050/prntbs.f \
	../m6050/prnvlm.f \
	../m6050/putone.f \
	../m6050/sumvw.f \
	../m6050/svlam.f \
	../m6050/swlam.f \
	../m6050/vlm.f \
	../m6050/wrtvlm.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6050

$(BINDIR)/m6050: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB) -o $(BINDIR)/m6050




