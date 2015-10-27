#   @(#)Makefile	1.1   11/30/90
#
# Makefile for M313
#
FC = f77
FFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../m313/dout32.f \
	../m313/dout64.f \
	../m313/dprim.f \
	../m313/driver.f \
	../m313/dsort.f \
	../m313/fixup.f \
	../m313/fmdint.f \
	../m313/fmint.f \
	../m313/getalp.f \
	../m313/hf2pdm.f \
	../m313/ldexp.f \
	../m313/m313.f \
	../m313/makdg.f \
	../m313/makdh.f \
	../m313/makdi.f \
	../m313/mkdghi.f \
	../m313/movder.f \
	../m313/noshel.f \
	../m313/nxtang.f \
	../m313/pikang.f \
	../m313/pikcen.f \
	../m313/pm313.f \
	../m313/prefac.f \
	../m313/prims.f \
	../m313/redund.f \
	../m313/vmakg.f \
	../m313/vmakh.f \
	../m313/vmaki2.f \
	../m313/vmkghi.f 


.f.o:
	$(FC) $(FFLAGS) -c $<

all: $(BINDIR)/m313

$(BINDIR)/m313: $(SRCS:.f=.o) $(MESALIB)
	$(FC) $(FFLAGS) $(SRCS:.f=.o) $(MESALIB) -o $(BINDIR)/m313

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m313.f
	cat $(SRCS) > ../source/m313.f

print: link
	lpr ../source/m313.f

clean:
	rm -f *.o core


