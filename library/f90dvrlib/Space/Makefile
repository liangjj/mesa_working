#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
DVRLIB = ../dvrlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
CPPFLAGS = -DDECPOINTER
FFLAGS = -c -O4 $(CPPFLAGS)
CO = co -l
RCS = RCS
REL = 

CPPSRCS = \
	blmat.F \
	cheb1.F \
	cheb2.F \
	conbsis.F \
	conham.F \
	dfadfb.F \
	diarep.F \
	dvrmul.F \
	faddfb.F \
	fadfb.F \
	fafb.F \
	filfun.F \
	filxv1.F \
	fn.F \
	hamphy.F \
	hermit.F \
	jacobi.F \
	kemat.F \
	lagure.F \
	lobatto.F \
	makd.F \
	maket.F \
	makod.F \
	ovmat.F \
	plymat.F \
	rmamp.F \
	tbloch.F \
	wtfn.F

.F.o:
	$(FC) $(FFLAGS) $<

all: $(DVRLIB) 
$(DVRLIB): $(CPPSRCS:.F=.o) 
	$(AR) $(ARFLAGS) $(DVRLIB) $(?:.F=.o)

sources: $(CPPSRCS)
	$CO) $(RCS)/*.F*