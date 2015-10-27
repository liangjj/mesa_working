#   %W%   %G%
#
#  Makefile for the LOCALBLAS routines
#
TOP_DIR = $(MESALIB)
include $(TOP_DIR)/Makefile.inc
.SUFFIXES: .f90 .mod $(SUFFIXES)

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
	../LOCALBLAS/dztpsv.f90 \
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
OBJECTS = \
        $(SRCS:.f=.o) \
        $(SRCS:.f90=.o)

all: $(BLASLIB)
$(BLASLIB): $(OBJECTS)
	$(AR) $(ARFLAGS) $(BLASLIB) $(OBJECTS)
	$(RANLIB) $(BLASLIB)
	$(AR) $(ARFLAGS) $(MESA_LIB) $(OBJECTS)
	$(RANLIB) $(MESA_LIB)

sources: $(SRCS)
$(SRCS):
	$(GET) $(RELEASE) $@
