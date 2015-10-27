#   %W%   %G%
#
#  Makefile for the LOCALBLAS routines
#
BLASLIB = ../../blaslib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c
GET = sccs get
REL = 

SRCS = \
	../LOCALBLAS/dasum.f \
	../LOCALBLAS/daxpy.f \
	../LOCALBLAS/dcopy.f \
	../LOCALBLAS/ddot.f \
	../LOCALBLAS/dgemm.f \
	../LOCALBLAS/dgemv.f \
	../LOCALBLAS/dscal.f \
	../LOCALBLAS/dswap.f \
	../LOCALBLAS/dsyr2k.f \
	../LOCALBLAS/dtrsm.f \
	../LOCALBLAS/idamax.f \
	../LOCALBLAS/izamax.f \
	../LOCALBLAS/lsame.f \
	../LOCALBLAS/zaxpy.f \
	../LOCALBLAS/zcopy.f \
	../LOCALBLAS/zdotc.f \
	../LOCALBLAS/zdotu.f \
	../LOCALBLAS/zgemm.f \
	../LOCALBLAS/zgemv.f \
	../LOCALBLAS/zgeru.f \
	../LOCALBLAS/zscal.f \
	../LOCALBLAS/zswap.f \
	../LOCALBLAS/ztrsm.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BLASLIB)

$(BLASLIB): $(SRCS:.f=.o) $(FOBJ)
	$(AR) $(ARFLAGS) $(BLASLIB) $(?:.f=.o)
	$(RANLIB) $(BLASLIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@
