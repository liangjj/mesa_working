#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the potential library
#
POTLIB = ../potlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

#CPPFLAGS = -DDECPOINTER
#FFLAGS = -c -O4 $(CPPFLAGS)
FFLAGS = -c -O4
REL = 

CPPSRCS = \
	addang.F \
	none.F \
	rwell.F \
	vadex12.F \
	vadpr.F \
	vanhmo.F \
	vcad12.F \
	vcoul.F \
	vexp.F \
	vexp12.F \
	vexp123.F \
	vhmo.F \
	vir4.F \
	vone.F \
	vpert.F \
	vres.F \
	vrwell.F \
	vwad12.F \
	vwell.F \
	vwel12.F \
	vwel123.F

.F.o:
	$(FC) $(FFLAGS) $<

all: $(POTLIB)
$(POTLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(POTLIB) $(?:.F=.o)

sources: $(CPPSRCS)
	co -l RCS/*.F
clean:
	rm *.o
