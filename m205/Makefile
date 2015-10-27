#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M205
#

GET = sccs get
REL = 

SRCS = \
	../m205/d2ecfrm.f \
	../m205/d2ecstp.f \
	../m205/d2eform.f \
	../m205/d2einfo.f \
	../m205/d2einit.f \
	../m205/d2emain.f \
	../m205/d2estep.f \
	../m205/eigint.f \
	../m205/fillfc.f \
	../m205/freq.f \
	../m205/m205.f \
	../m205/pm205.f \
	../m205/tvec.f \
	../m205/vibfrc.f \
	../m205/vibfreq.f \
	../m205/zctran.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m205

$(BINDIR)/m205: $(SRCS:.f=.o) $(MESALIB)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(MESALIB) $(MDLIB) -o $(BINDIR)/m205

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m205.f
	cat $(SRCS) > ../source/m205.f

print: link
	lpr ../source/m205.f

clean:
	rm -f *.o core



