#   %W%   %G%
#
# Makefile for uass
#
# FC = g77
FFLAGS = $(FFLAGSOZ)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
BLAS= $(BLASUSE)
GET = sccs get
REL = 

SRCS = \
	../kuass/lpmn.f

SRCSV = \
	../kuass/uass2.for 

all: $(BINDIR)/xass

$(SRCS:%.f=%.o): %.o: %.f
	$(FC) $(FFLAGS) $<

$(SRCSV:%.for=%.o): %.o: %.for
	sed -e "s/#maxnbfkohn/$(MAXNBFKOHN)/g" <$(*).for | \
	sed -e "s/#maxsmall/$(MAXSMALL)/g" | \
	sed -e "s/#maxbig/$(MAXBIG)/g" | \
	sed -e "s/#maxltop/$(MAXLTOP)/g" | \
	sed -e "s/#maxchan/$(MAXCHAN)/g" >$(*).f
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xass:  $(SRCSV:%.for=%.o) $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS)  $(SRCSV:%.for=%.o) $(SRCS:.f=.o) $(MESALIB) $(BLAS) -o $(BINDIR)/xass

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/uass2.f
	cat $(SRCS) > ../source/uass2.f

print: link
	lpr ../source/uass2.f

clean:
	rm -f *.o core
