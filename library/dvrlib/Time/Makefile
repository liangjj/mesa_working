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
CO = -co -l
RCS = RCS
REL = 

CPPSRCS = \
	fdf.F \
	nrmlze.F \
	timpts.F \
	tlobato.F \
	tphy.F \
	tplmat.F \
	twpt.F

.F.o:
	$(FC) $(FFLAGS) $<

all: $(DVRLIB)
$(DVRLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(DVRLIB) $(?:.F=.o)

sources: $(CPPSRCS)
	$(CO) $(RCS)/*.F*

clean:
	rm *.o

