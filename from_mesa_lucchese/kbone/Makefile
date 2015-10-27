#   %W%   %G%
#
# Makefile for bone
#
# FC = g77
FFLAGS = $(FFLAGSZ)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
GET = sccs get
REL = 

SRCS = \
	../kbone/class.f \
	../kbone/cubfev.f \
	../kbone/gamfun.f \
	../kbone/gaussq.f \
	../kbone/gbslve.f \
	../kbone/gbtql2.f \
	../kbone/monder.f \
	../kbone/plm.f \
	../kbone/pts.f \
	../kbone/pwcfev.f \
	../kbone/radial.f \
	../kbone/sjymec.f 

SRCSV = \
	../kbone/bnhead.for

all: $(BINDIR)/xbone

$(SRCS:%.f=%.o): %.o: %.f
	$(FC) $(FFLAGS) $<

$(SRCSV:%.for=%.o): %.o: %.for
	sed -e "s/#maxltop/$(MAXLTOP)/g" <$(*).for >$(*).f
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xbone: $(SRCSV:%.for=%.o) $(SRCS:.f=.o) $(MESALIB) $(BLASLIB)
	$(LD) $(LDFLAGS) $(SRCSV:%.for=%.o) $(SRCS:.f=.o)  -o $(BINDIR)/xbone

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/bone.f
	cat $(SRCS) > ../source/bone.f

print: link
	lpr ../source/bone.f

clean:
	rm -f *.o core
