#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M7500
#
FC = f77
FFLAGS = -c -O4

LD = f77
LDFLAGS =
MAIN = /usr1/mesa
#MAIN = /usr1/users/bis/mesa
BINDIR = $(MAIN)/bin
MESCTLIB = $(MAIN)/library/mesctlib.a
VECLIB = $(MAIN)/library/libdxml.a

GET = sccs get

SRCS = \
	gauwpt.f \
	genq.f \
	genrwt.f \
	lancz.f \
	matprnt.f \
	orpoly.f \
	scaprd.f \
	sumary.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m7500

$(BINDIR)/m7500: $(SRCS:.f=.o) $(MESCTLIB) $(VECLIB)  
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(MESCTLIB) $(VECLIB) \
					-o $(BINDIR)/m7500

