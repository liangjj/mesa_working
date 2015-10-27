#   %W%   %G%
#
# Makefile for mfpad
#
# FC = g77
FFLAGS = $(FFLAGSOZ)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
#MESALIB = ../library/mesalib.a
#BLASLIB= $(BLASUSE)
GET = sccs get
REL = 

SRCS = \
	../kmfpad/mfpad.for

all: $(BINDIR)/xmfpad

$(SRCS:%.for=%.o): %.o: %.for
	sed -e "s/#maxltop/$(MAXLTOP)/g" <$(*).for >$(*).f
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xmfpad: $(SRCS:%.for=%.o) 
	$(LD) $(LDFLAGS) $(SRCS:%.for=%.o) -o $(BINDIR)/xmfpad

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/mfpad.f
	cat $(SRCS) > ../source/mfpad.f

print: link
	lpr ../source/mfpad.f

clean:
	rm -f *.o core
