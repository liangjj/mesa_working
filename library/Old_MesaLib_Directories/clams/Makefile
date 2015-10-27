
#   %W%   %G%
#
#  Makefile for the CLAMS (common Los Alamos mathematical subroutines)
#  portion of MesaLib
#
MYLIB = ../mylib.a
MESALIB = ../mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

REL = 

CPPSRCS = \
	../clams/caxpy.F \
	../clams/ccopy.F \
	../clams/cdotc.F \
	../clams/cdotu.F \
	../clams/cgefa.F \
	../clams/cgemm.F \
	../clams/cgemv.F \
	../clams/cgesl.F \
	../clams/cgetrf.F \
	../clams/cgetrs.F \
	../clams/cscal.F \
	../clams/cswap.F \
	../clams/epslon.F \
	../clams/icamax.F \
	../clams/icopy.F \
	../clams/imtql1.F \
	../clams/imtql2.F \
	../clams/isamax.F \
	../clams/ismin.F \
	../clams/mxma.F \
	../clams/pythag.F \
	../clams/rsp.F \
	../clams/sasum.F \
	../clams/saxpy.F \
	../clams/scatter.F \
	../clams/scopy.F \
	../clams/sdot.F \
	../clams/sgeco.F \
	../clams/sgedi.F \
	../clams/sgefa.F \
	../clams/sgefs.F \
	../clams/sgemm.F \
	../clams/sgemv.F \
	../clams/sgesl.F \
	../clams/sgetrf.F \
	../clams/sgetrs.F \
	../clams/sscal.F \
	../clams/sspfa.F \
	../clams/sspsl.F \
	../clams/sswap.F \
	../clams/tql1.F \
	../clams/tql2.F \
	../clams/tqlrat.F \
	../clams/trbak1.F \
	../clams/trbak3.F \
	../clams/tred1.F \
	../clams/tred2.F \
	../clams/tred3.F 

.F.o:
	$(FC) $(FFLAGS) $(CPPFLAGS) $<

all: $(MYLIB)
$(MYLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(MYLIB) $(?:.F=.o)
	$(RANLIB) $(MYLIB)
all: $(MESALIB)
$(MESALIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(MESALIB) $(?:.F=.o)
	$(RANLIB) $(MESALIB)

sources: $(CPPSRCS)
$(CPPSRCS):
	co -l RCS/*.F*

clean: $(CPPSRCS)
	rm *.F *.o



