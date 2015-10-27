#   %W%   %G%
#
# Makefile for mfpad
#
FC = g77
FFLAGS = -c -O 

LD = g77
LDFLAGS =

BINDIR = ../bin
#BLASLIB= /System/Library/Frameworks/vecLib.framework/versions/A/vecLib
REL = 

SRCS = mfpad3D.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/xmfpad3D

$(BINDIR)/xmfpad3D: $(SRCS:.f=.o) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) -o $(BINDIR)/xmfpad3D



clean:
	rm -f *.o core $(BINDIR)/xmfpad3D 
