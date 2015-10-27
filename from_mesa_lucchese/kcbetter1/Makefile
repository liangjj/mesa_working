#   %W%   %G%
#
# Makefile for cbetter1
#
# FC = g77
FFLAGS = $(FFLAGSO)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
#BLASLIB= $(BLASUSE)
GET = sccs get
REL = 

SRCS = \
	../kcbetter1/cbetter1.for

all: $(BINDIR)/xcbetter1

$(SRCS:%.for=%.o): %.o: %.for
	sed -e "s/#maxnbfkohn/$(MAXNBFKOHN)/g" <$(*).for | \
	sed -e "s/#maxltop/$(MAXLTOP)/g" | \
	sed -e "s/#maxchan/$(MAXCHAN)/g" >$(*).f
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xcbetter1: $(SRCS:%.for=%.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:%.for=%.o) $(MESALIB) -o $(BINDIR)/xcbetter1

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/cbetter1.f
	cat $(SRCS) > ../source/cbetter1.f

print: link
	lpr ../source/cbetter1.f

clean:
	rm -f *.o core
