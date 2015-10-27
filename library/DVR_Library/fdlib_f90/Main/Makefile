#   @(#)Makefile	1.4   8/2/91
#
# Makefile for TIMLIB
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR=$(DVR)
include $(TOP_DIR)/Makefile.inc
FD_LIB = ../../fd_lib_f90.a
RCS = RCS
CO = co -l

SRCS = \
	band3.f90 \
        band5.f90 \
	band7.f90 \
	cpy_3.f90 \
	fd_input.f90 \
	fd_basis.f90 \
	points.f90 


all: $(FD_LIB)
$(FD_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(FD_LIB) $(?:.f90=.o)

sources: $(SRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o
