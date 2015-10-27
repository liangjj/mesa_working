#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M319
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
	../m319/cfrt1.f \
	../m319/driver.f \
	../m319/m319.f \
	../m319/pm319.f 


.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m319

$(BINDIR)/m319: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m319

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m319.f
	cat $(SRCS) > ../source/m319.f

print: link
	lpr ../source/m319.f

clean:
	rm -f *.o core


