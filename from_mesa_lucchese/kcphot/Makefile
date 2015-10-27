#   %W%   %G%
#
# Makefile for cphot
#
# FC = g77
FFLAGS = $(FFLAGSOZ)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
#BLASLIB= $(BLASUSE)
GET = sccs get
REL = 

SRCSV = \
	../kcphot/cphot.for

all: $(BINDIR)/xcphot

$(SRCSV:%.for=%.o): %.o: %.for
	sed -e "s/#maxnbfkohn/$(MAXNBFKOHN)/g" <$(*).for | \
	sed -e "s/#maxsmall/$(MAXSMALL)/g" | \
	sed -e "s/#maxlmtop/$(MAXLMTOP)/g" | \
	sed -e "s/#maxchan/$(MAXCHAN)/g" >$(*).f
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xcphot: $(SRCSV:%.for=%.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCSV:%.for=%.o) $(MESALIB) -o $(BINDIR)/xcphot


clean:
	rm -f *.o core
