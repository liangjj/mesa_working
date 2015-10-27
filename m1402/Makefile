#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M1402
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = $(MESA_BIN)
MESALIB = $(MESA_LIB)/mesalib.a
MDLIB =$(MDLIB)
GET = sccs get
REL = 

SRCS = \
	../m1402/fixden.f \
	../m1402/m1402.f \
	../m1402/pm1402.f \
	../m1402/vaeapb.f \
	../m1402/vecout.f \
	../m1402/vseamb.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m1402

$(BINDIR)/m1402: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m1402

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m1402.f
	cat $(SRCS) > ../source/m1402.f

print: link
	lpr ../source/m1402.f

clean:
	rm -f *.o core


