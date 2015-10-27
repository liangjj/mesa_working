#   @(#)Makefile	5.1   11/6/94
#
#
FC = ifort
LD = $(FC)
LDFLAGS =

BINDIR = $(MESA_BIN)
MESALIB = $(MESA_LIB)/mesalib.a
GET = sccs get
REL = 

SRCS = \
	../m999/test_type.f90 \
	../m101/cprint.f \
	../m101/doc101.f \
	../m101/isubst.f \
	../m101/lsubst.f \
	../m101/m101.f \
	../m101/pm101.f \
	../m101/zcentr.f \
	../m101/zget.f \
	../m101/zmatch.f \
	../m101/zparm.f \
	../m101/zsymbl.f \
	../m101/zvar.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m101

$(BINDIR)/m101: $(SRCS:.f=.o) $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m101

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m101.f
	cat $(SRCS) > ../source/m101.f

print: link
	lpr ../source/m101.f

clean:
	rm -f *.o core
