#  @(#)Makefile	1.8 8/1/91
#
# Makefile for M6202
#
FC = f77
FFLAGS = -c +OP4  

LD = f77
LDFLAGS =

BINDIR = ../bin
MESCTLIB = ../library/mesctlib.a
VECLIB = /usr/lib/libvec.a
GET = sccs get
REL = 

SRCS = \
	caspln.f \
	dgsn.f \
	ddgsn.f \
	ftoylm.f \
	gsn.f \
	legend.f \
	legint.f \
	lmdcmp.f \
	miscfn.f \
	mkspln.f \
	mkylm.f \
	nrmylm.f \
	plmtof.f \
	runiqe.f \
	scmphi.f \
	setzro.f \
	sqrtwt.f \
	srtleg.f \
	ylmint.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6202

$(BINDIR)/m6202: $(SRCS:.f=.o)  ${MESCTLIB} $(VECLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) ${MESCTLIB} $(VECLIB) -o $(BINDIR)/m6202
