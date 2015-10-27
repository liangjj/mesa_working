#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M994
#
FC = fort 
FFLAGS = -c 
LD = fort
LDFLAGS =

BINDIR = ../m994
GET = sccs get
REL = 

SRCS = \
	m994.f \
	pm994.f \
	iostd.f \
	readran.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m994

$(BINDIR)/m994: $(SRCS:.f=.o) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) -o $(BINDIR)/m994
