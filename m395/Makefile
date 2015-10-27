#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M395
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
	../m395/drvint.f \
	../m395/gtrans.f \
	../m395/m395.f \
	../m395/prod.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m395

$(BINDIR)/m395: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m395

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m395.f
	cat $(SRCS) > ../source/m395.f

print: link
	lpr ../source/m395.f

clean:
	rm -f *.o core



