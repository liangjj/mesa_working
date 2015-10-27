#   @(#)Makefile	1.9   8/7/91
#
#
FC = f77
FFLAGS = -c -O4
LD = f77
LDFLAGS =
MAIN = /usr1/mesa
BINDIR = $(MAIN)/bin
MESCTLIB = $(MAIN)/library/mesctlib.a
COULIB = $(MAIN)/library/coulib/coulib.a

GET = sccs get
REL = 

SRCS = \
	asycou.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m6297

$(BINDIR)/m6297: $(SRCS:.f=.o) $(COULIB) $(MESCTLIB)   
	$(LD) $(LDFLAGS) $(SRCS:.f=.o)  $(COULIB) $(MESCTLIB) \
					-o $(BINDIR)/m6297


