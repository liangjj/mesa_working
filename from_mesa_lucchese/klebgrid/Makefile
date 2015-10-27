#   %W%   %G%
#
# Makefile for Lebgrid
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
	../klebgrid/class.f \
	../klebgrid/gamfun.f \
	../klebgrid/gaussq.f \
	../klebgrid/gbslve.f \
	../klebgrid/gbtql2.f \
	../klebgrid/lebgrid.f \
	../klebgrid/voronoi.f 

SRCSF90 = \
	../klebgrid/LebedevLaikov.f90

all: $(BINDIR)/xlebgrid

$(SRCS:%.f=%.o): %.o: %.f
	$(FC) $(FFLAGS) $<

$(SRCSF90:%.f90=%.o): %.o: %.f90
	$(FC90) $(FFLAGS90) -c $(*).f90

$(BINDIR)/xlebgrid: $(SRCS:.f=.o) $(SRCSF90:%.f90=%.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(SRCSF90:%.f90=%.o) $(MESALIB) -o $(BINDIR)/xlebgrid

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/lebgrid.f
	cat $(SRCS) > ../source/lebgrid.f

print: link
	lpr ../source/lebgrid.f

clean:
	rm -f *.o core
