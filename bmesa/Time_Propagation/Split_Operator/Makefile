#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M8015
#
FC = f77
FFLAGS = -c -O4

LD = f77
LDFLAGS =
MAIN = /usr1/mesa
SUB = $(MAIN)/bmesa/Time_Propagation
LIBDIR = $(MAIN)/library
BINDIR = $(MAIN)/bin
DVRLIB = $(LIBDIR)/dvrlib/dvrlib.a
POTLIB = $(LIBDIR)/potlib/potlib.a
HAMLIB = $(LIBDIR)/hamlib/hamlib.a
MESCTLIB = $(LIBDIR)/mesctlib.a
ITDVRLIB  = $(LIBDIR)/itlib/itdvrlib.a
TDVRLIB = $(LIBDIR)/tdvrlib/tdvrlib.a
TIMLIB = $(SUB)/timlib.a
VECLIB = $(LIBDIR)/libdxml.a

GET = sccs get

SRCS = \
	cvexp.f \
	drivso.f \
	ke1exp.f \
	ke2exp.f \
	ke3exp.f \
	prp3d.f \
	so3d.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m8015

$(BINDIR)/m8015: $(SRCS:.f=.o)	$(TIMLIB) $(DVRLIB) $(HAMLIB) $(POTLIB) \
				$(ITDVRLIB) $(TDVRLIB) $(MESCTLIB) $(VECLIB)

				$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(TIMLIB) $(DVRLIB) \
				$(HAMLIB) $(POTLIB) $(ITDVRLIB) $(TDVRLIB) \
				$(MESCTLIB) $(VECLIB) \
				-o $(BINDIR)/m8015
