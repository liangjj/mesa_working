#   %W%   %G%
#
# Makefile for psff
#
# FC = g77
FFLAGS = $(FFLAGSOZ)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MYLIB = ../dir_mylib/mylib.a
BLASLIB= $(BLASUSE)
GET = sccs get
REL = 

SRCS = \
	../kpsff/cmwrite.f \
	../kpsff/fun.f \
	../kpsff/matinvr.f \
	../kpsff/plm.f \
	../kpsff/qagi.f \
	../kpsff/qagie.f \
	../kpsff/qelg.f \
	../kpsff/qk15i.f \
	../kpsff/qpsrt.f \
	../kpsff/spline.f 

SRCSV = \
	../kpsff/psff.for

all: $(BINDIR)/xpsff

$(SRCS:%.f=%.o): %.o: %.f
	$(FC) $(FFLAGS) $<

$(SRCSV:%.for=%.o): %.o: %.for
	sed -e "s/#maxnbfkohn/$(MAXNBFKOHN)/g" <$(*).for | \
	sed -e "s/#maxltop/$(MAXLTOP)/g" | \
	sed -e "s/#maxlmtop/$(MAXLMTOP)/g" | \
	sed -e "s/#maxchan/$(MAXCHAN)/g" >$(*).f
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xpsff: $(SRCSV:%.for=%.o) $(SRCS:.f=.o)
	$(LD) $(LDFLAGS) $(SRCSV:%.for=%.o) $(SRCS:.f=.o) $(MYLIB) $(MESALIB) $(BLASLIB) -o $(BINDIR)/xpsff

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/psff.f
	cat $(SRCS) > ../source/psff.f

print: link
	lpr ../source/psff.f

clean:
	rm -f *.o core
