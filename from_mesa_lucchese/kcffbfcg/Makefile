#   %W%   %G%
#
# Makefile for cffbfcg
#
# FC = g77
FFLAGS = $(FFLAGSOZ)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
BLASLIB  =$(BLASUSE)
GET = sccs get
REL = 

SRCS = \
	../kcffbfcg/bfints.f \
	../kcffbfcg/ffints.f 

SRCSV = \
	../kcffbfcg/cffbfcg.for

all: $(BINDIR)/xcffbfcg

$(SRCS:%.f=%.o): %.o: %.f
	$(FC) $(FFLAGS) $<

$(SRCSV:%.for=%.o): %.o: %.for
	sed -e "s/#maxnbfkohn/$(MAXNBFKOHN)/g" <$(*).for | \
	sed -e "s/#maxltop/$(MAXLTOP)/g" | \
	sed -e "s/#maxlmtop/$(MAXLMTOP)/g" | \
	sed -e "s/#maxchan/$(MAXCHAN)/g" >$(*).f
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xcffbfcg: $(SRCSV:%.for=%.o) $(SRCS:.f=.o)
	$(LD) $(LDFLAGS) $(SRCSV:%.for=%.o) $(SRCS:.f=.o) $(MESALIB) $(BLASLIB) -o $(BINDIR)/xcffbfcg

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/cffbfcg.f
	cat $(SRCS) > ../source/cffbfcg.f

print: link
	lpr ../source/cffbfcg.f

clean:
	rm -f *.o core
