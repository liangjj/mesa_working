#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M8011
#
FC = f90
FFLAGS = -c -O4

LD = f90
LDFLAGS =
MAIN = /usr1/mesa
LIBDIR = $(MAIN)/library
BINDIR = ./
LAPACK  = $(LIBDIR)/LAPACK/lapack.a
LIBDXML = $(LIBDIR)/libdxml.a

GET = sccs get

SRCS = \
	spop_1d.f90

.f90.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m8011

$(BINDIR)/m8011: $(SRCS:.f90=.o) $(LAPACK) $(LIBDXML)
	$(LD) $(LDFLAGS) $(SRCS:.f90=.o)  $(LAPACK) $(LIBDXML) \
					-o $(BINDIR)/m8011




