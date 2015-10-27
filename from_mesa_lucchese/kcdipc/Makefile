#   %W%   %G%
#
# Makefile for cdipc
#
# FC = g77
FFLAGS = $(FFLAGSOZ)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
BLASLIB= $(BLASUSE)
GET = sccs get
REL = 

SRCS = \
     ../kcdipc/cdipint.f \
     ../kcdipc/bforthog.f \
     ../kcdipc/trnmat.f \
     ../kcdipc/rdbinsqr.f \
     ../kcdipc/rdbintri.f

SRCSV = \
	../kcdipc/cdipc.for

all: $(BINDIR)/xcdipc

$(SRCS:%.f=%.o): %.o: %.f
	$(FC) $(FFLAGS) $<

$(SRCSV:%.for=%.o): %.o: %.for
	sed -e "s/#maxnbfkohn/$(MAXNBFKOHN)/g" <$(*).for | \
	sed -e "s/#maxltop/$(MAXLTOP)/g" | \
	sed -e "s/#maxlmtop/$(MAXLMTOP)/g" | \
	sed -e "s/#maxchan/$(MAXCHAN)/g" >$(*).f
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xcdipc: $(SRCSV:%.for=%.o) $(SRCS:.f=.o)
	$(LD) $(LDFLAGS) $(SRCSV:%.for=%.o) $(SRCS:.f=.o) $(MESALIB) $(BLASLIB) -o $(BINDIR)/xcdipc

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/cdipc.f
	cat $(SRCS) > ../source/cdipc.f

print: link
	lpr ../source/cdipc.f

clean:
	rm -f *.o core
