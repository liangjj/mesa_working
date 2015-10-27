#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the dvr routines
#
.SUFFIXES: .f90 $(SUFFIXES)
DVR_LIB = ../../dvr_lib_f90.a
TOP_DIR=$(DVR)
include $(TOP_DIR)/Makefile.inc
RCS = RCS
CO = co -l
REL = 

SRCS = \
	blcyl.f90 \
	blrad.f90 \
	blxyz.f90 \
	cnvtab.f90 \
	cnvtpt.f90 \
	genwtf.f90 \
	genq.f90 \
	genrwt.f90 \
	kemel.f90 \
	kinxyz.f90 \
	kinrad.f90 \
	kincyl.f90 \
	lancz.f90 \
	modab.f90 \
	plywpt.f90 \
	recur.f90 \
	scaprd.f90 \
	todelta.f90 \
	xsq2x.f90

all: $(DVR_LIB)
$(DVR_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(DVR_LIB) $(?:.f90=.o)
	$(RANLIB) $(DVR_LIB)
sources: $(SRCS)
	$(CO) $(RCS)/*.f90*

clean:
	rm -f *.o

