#   @(#)Makefile	1.4   8/2/91
#
# Makefile for TPRPLIB
#
.SUFFIXES: .f90 .mod $(SUFFIXES)
RANLIB = ranlib
TIMPRP = $(MESA_HOME)/bmesa/Time_Propagation
TPRPLIB = $(TIMPRP)/tprplib_f90.a

GET = sccs get

SRCS = \
	lanczos_global.f90 \
	lanczos_prnt.f90 \
	lanczos_dat.f90 \
	lanczos.f90 \
	lanczos_prop.f90 \
	tri_diag.f90

.f90.o:
	$(FC) $(FFLAGS) -I$(MOD_DIR) $<

all: $(TPRPLIB)
$(TPRPLIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(TPRPLIB) $(?:.f90=.o)
