#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M992
#
FC = fort 
FFLAGS = -c 
LD = fort
LDFLAGS =

BINDIR = ../m992
GET = sccs get
REL = 

SRCS = \
	m992.f \
	pm992.f \
	iostd.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m992

$(BINDIR)/m992: $(SRCS:.f=.o) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) -o $(BINDIR)/m992

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m992.f
	cat $(SRCS) > ../source/m992.f

print: link
	lpr ../source/m992.f

clean:
	rm -f *.o core



