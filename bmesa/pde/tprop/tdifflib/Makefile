#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the differential equation routines
#
TDIFFLIB = tdifflib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = f77 
FFLAGS = -c 
LD = f77
LDFLAGS =
GET = sccs get
REL = 

SRCS = \
	deabm.f \
	des.f \
	dstrt.f \
	dhvnm.f \
	dntp.f \
	dstps.f \
	odeint.f \
	rk.f \
	rk4.f \
	rkqs.f \
	rkck.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(TDIFFLIB)
$(TDIFFLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(TDIFFLIB) $(?:.f=.o)


