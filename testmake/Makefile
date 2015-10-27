#   @(#)Makefile	1.4   8/1/91
#
# Makefile for M940
#
AR = -ar
ARFLAGS = rvu
FC = fort
CPPFLAGS = -DDECFLAG
FFLAGS = -c -O5 $(CPPFLAGS)
BINDIR = testmake

LD = fort
LDFLAGS =

TESTLIB = testlib.a
REL = 

SRCS = test.f

CPPSRCS = \
	test1.F

.F.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)
lib: $(TESTLIB)
$(BINDIR): $(SRCS:.f=.o) $(TESTLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(TESTLIB) -o $(BINDIR)
$(TESTLIB): $(CPPSRCS:.F=.o)
	$(AR) $(ARFLAGS) $(TESTLIB) $(?:.F=.o)

sources: $(SRCS) $(CPPSRCS)
$(SRCS):
	co -l RCS/*.f*
$(CPPSRCS):
	co -l RCS/*.F*



