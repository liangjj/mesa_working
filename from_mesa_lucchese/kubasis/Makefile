#   %W%   %G%
#
# Makefile for ubasis
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
	../kubasis/gaussval.f \
	../kubasis/inputs.f

SRCSV = \
	../kubasis/ubasis.for 

all: $(BINDIR)/xbasis

$(SRCS:%.f=%.o): %.o: %.f
	$(FC) $(FFLAGS) $<

$(SRCSV:%.for=%.o): %.o: %.for
	sed -e "s/#maxnbfkohn/$(MAXNBFKOHN)/g" <$(*).for | \
	sed -e "s/#maxchan/$(MAXCHAN)/g" | \
	sed -e "s/#maxprimkohn/$(MAXPRIMKOHN)/g" >$(*).f	  
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xbasis: $(SRCSV:%.for=%.o) $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCSV:%.for=%.o) $(SRCS:.f=.o) $(MESALIB) $(BLAS) -o $(BINDIR)/xbasis

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/ubasis.f
	cat $(SRCS) > ../source/ubasis.f

print: link
	lpr ../source/ubasis.f

clean:
	rm -f *.o core
