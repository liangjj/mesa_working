#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M323
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS = 

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB =
GET = sccs get
REL = 

SRCS = \
	../m323/dout32.f \
	../m323/dout64.f \
	../m323/dprim.f \
	../m323/driver.f \
	../m323/dsort.f \
	../m323/fixup.f \
	../m323/fmdint.f \
	../m323/fmint.f \
	../m323/getalp.f \
	../m323/ldexp.f \
	../m323/m323.f \
	../m323/makdg.f \
	../m323/makdh.f \
	../m323/makdi.f \
	../m323/mkdghi.f \
	../m323/movder.f \
	../m323/noshel.f \
	../m323/nxtang.f \
	../m323/pikang.f \
	../m323/pikcen.f \
	../m323/pm323.f \
	../m323/prefac.f \
	../m323/prims.f \
	../m323/redund.f \
	../m323/vmakg.f \
	../m323/vmakh.f \
	../m323/vmaki2.f \
	../m323/vmkghi.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m323

$(BINDIR)/m323: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m323

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m323.f
	cat $(SRCS) > ../source/m323.f

print: link
	lpr ../source/m323.f

clean:
	rm -f *.o core


