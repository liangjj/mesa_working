#   %W%   %G%
#
# Makefile for ustatnuc
#
# FC = g77
FFLAGS = $(FFLAGSZ)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
BLAS= $(BLASUSE)
GET = sccs get
REL = 

SRCS = \
	../kustatnuc/elim.f \
	../kustatnuc/fmch.f \
	../kustatnuc/generf.f \
	../kustatnuc/gfunct.f \
	../kustatnuc/inputs.f \
	../kustatnuc/stuff0.f \
	../kustatnuc/stuff1.f \
	../kustatnuc/stuff2.f \
	../kustatnuc/stuff3.f \
	../kustatnuc/stuff4.f \
	../kustatnuc/stuff5.f \
	../kustatnuc/stuff6.f \
	../kustatnuc/rho.f \
	../kustatnuc/vints.f 

SRCSV = \
	../kustatnuc/ustatnuc.for

all: $(BINDIR)/xstatnuc

$(SRCS:%.f=%.o): %.o: %.f
	$(FC) $(FFLAGS) $<

$(SRCSV:%.for=%.o): %.o: %.for
	sed -e "s/#maxnbfkohn/$(MAXNBFKOHN)/g" <$(*).for | \
	sed -e "s/#maxchan/$(MAXCHAN)/g" | \
	sed -e "s/#maxprimkohn/$(MAXPRIMKOHN)/g" >$(*).f	  
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xstatnuc: $(SRCSV:%.for=%.o) $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCSV:%.for=%.o) $(SRCS:.f=.o) $(MESALIB) $(BLAS) -o $(BINDIR)/xstatnuc

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f print.f
	cat $(SRCS) > print.f

print: link
	lpr print.f

clean:
	rm -f *.o core

