#   %W%   %G%
#
# Makefile for csolve
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
	../kcsolve/cgeco.f \
	../kcsolve/cgefa.f \
	../kcsolve/cgefs.f \
	../kcsolve/cgesl.f \
	../kcsolve/expand.f \
	../kcsolve/fdump.f \
	../kcsolve/i1mach.f \
	../kcsolve/j4save.f \
	../kcsolve/r1mach.f \
	../kcsolve/rmachr4.f \
	../kcsolve/xerprn.f \
	../kcsolve/xerhlt.f \
	../kcsolve/xercnt.f \
	../kcsolve/xersve.f \
	../kcsolve/xgetua.f \
	../kcsolve/rdbinsqr.f \
	../kcsolve/xermsg.f 

SRCSV = \
	../kcsolve/csolve.for

all: $(BINDIR)/xcsolve

$(SRCS:%.f=%.o): %.o: %.f
	$(FC) $(FFLAGS) $<

$(SRCSV:%.for=%.o): %.o: %.for
	sed -e "s/#maxnbfkohn/$(MAXNBFKOHN)/g" <$(*).for | \
	sed -e "s/#maxsmall/$(MAXSMALL)/g" | \
	sed -e "s/#maxbig/$(MAXBIG)/g" | \
	sed -e "s/#maxchan/$(MAXCHAN)/g" >$(*).f
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xcsolve: $(SRCS:.f=.o) $(SRCSV:%.for=%.o)
	$(LD) $(LDFLAGS) $(SRCSV:%.for=%.o) $(SRCS:.f=.o) $(MESALIB) $(BLASLIB) -o $(BINDIR)/xcsolve

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/csolve.f
	cat $(SRCS) > ../source/csolve.f

print: link
	lpr ../source/csolve.f

clean:
	rm -f *.o core
