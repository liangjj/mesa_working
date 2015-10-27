#   %W%   %G%
#
# Makefile for mfpad
#
# FC = g77
FFLAGS = $(FFLAGSO)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
#BLASLIB= $(BLASUSE)
REL = 

SRCS = mfpad3D.for

all: $(BINDIR)/xmfpad3D

$(SRCS:%.for=%.o): %.o: %.for
	sed -e "s/#maxltop/$(MAXLTOP)/g" <$(*).for >$(*).f
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xmfpad3D: $(SRCS:%.for=%.o) 
	$(LD) $(LDFLAGS) $(SRCS:%.for=%.o) -o $(BINDIR)/xmfpad3D



clean:
	rm -f *.o core $(BINDIR)/xmfpad3D 
