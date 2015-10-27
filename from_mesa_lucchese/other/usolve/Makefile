#   %W%   %G%
#
# Makefile for usolve
#
FC = g77
FFLAGS = -c -O -finit-local-zero -finit-local-zero

LD = g77
LDFLAGS =

BINDIR = ../bin
BLASLIB= /System/Library/Frameworks/vecLib.framework/versions/A/vecLib
GET = sccs get
REL = 

SRCS = \
	../usolve/cgeco.f \
	../usolve/cgefa.f \
	../usolve/cgefs.f \
	../usolve/cgesl.f \
	../usolve/fdump.f \
	../usolve/i1mach.f \
	../usolve/j4save.f \
	../usolve/r1mach.f \
	../usolve/xerprn.f \
	../usolve/xerhlt.f \
	../usolve/xercnt.f \
	../usolve/xersve.f \
	../usolve/xgetua.f \
	../usolve/usolve.f \
	../usolve/rdbinsqr.f \
	../usolve/xermsg.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/xusolve

$(BINDIR)/xusolve: $(SRCS:.f=.o) $(BLASLIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(BLASLIB) -o $(BINDIR)/xusolve

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/usolve.f
	cat $(SRCS) > ../source/usolve.f

print: link
	lpr ../source/usolve.f

clean:
	rm -f *.o core
