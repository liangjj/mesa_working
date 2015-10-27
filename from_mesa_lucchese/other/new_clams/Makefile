#   %W%   %G%
#
#  Makefile for the CLAMS (common Los Alamos mathematical subroutines)
#  portion of MesaLib
#
MESALIB = ../mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = g77
FFLAGS = -c -O2
GET = sccs get
REL = 

SRCS = \
	../clams/imtql2.f \
	../clams/isamax.f \
	../clams/ismin.f \
	../clams/mxma.f \
	../clams/pythag.f \
	../clams/rsp.f \
	../clams/sasum.f \
	../clams/saxpy.f \
	../clams/scatter.f \
	../clams/scopy.f \
	../clams/sdot.f \
	../clams/sgeco.f \
	../clams/sgedi.f \
	../clams/sgefa.f \
	../clams/sgemm.f \
	../clams/sgemv.f \
	../clams/sgesl.f \
	../clams/sscal.f \
	../clams/sspev.f \
	../clams/sspfa.f \
	../clams/sspsl.f \
	../clams/sswap.f \
	../clams/tql2.f \
	../clams/tqlrat.f \
	../clams/trbak3.f \
	../clams/tred3.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(MESALIB)
$(MESALIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(MESALIB) $(?:.f=.o)
	$(RANLIB) $(MESALIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@

link: sources
	rm -f ../../source/clams.f
	cat $(SRCS) > ../../source/clams.f

print: link
	lpr ../../source/clams.f

clean:
	rm -f *.o
	rm -f core

