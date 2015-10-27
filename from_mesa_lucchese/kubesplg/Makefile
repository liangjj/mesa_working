#   %W%   %G%
#
# Makefile for ubesplg
#
# FC = g77
FFLAGS = $(FFLAGSOZ)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../kubesplg/outward.f \
	../kubesplg/sjymec.f \
	../kubesplg/spline.f \
	../kubesplg/splinec.f \
	../kubesplg/zero.f 

SRCSV = \
	../kubesplg/gg.for \
	../kubesplg/ubesplg.for

all: $(BINDIR)/xbesplg

$(SRCS:%.f=%.o): %.o: %.f
	$(FC) $(FFLAGS) $<

$(SRCSV:%.for=%.o): %.o: %.for
	sed -e "s/#maxltop/$(MAXLTOP)/g" <$(*).for | \
	sed -e "s/#maxchan/$(MAXCHAN)/g" >$(*).f
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xbesplg: $(SRCS:%.f=%.o) $(SRCSV:%.for=%.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:%.f=%.o) $(SRCSV:%.for=%.o) $(MESALIB) -o $(BINDIR)/xbesplg

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/ubesplg.f
	cat $(SRCS) > ../source/ubesplg.f

print: link
	lpr ../source/ubesplg.f

clean:
	rm -f *.o core
