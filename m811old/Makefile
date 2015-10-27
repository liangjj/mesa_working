#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M811
#
include $(MESA_HOME)/Makefile.inc
SOURCE = m811
GET = sccs get
REL = 

SRCS = \
	$(MESA_HOME)/$(SOURCE)old/chknrm.f \
	$(MESA_HOME)/$(SOURCE)old/fock.f \
	$(MESA_HOME)/$(SOURCE)old/getdrt.f \
	$(MESA_HOME)/$(SOURCE)old/$(SOURCE).f \
	$(MESA_HOME)/$(SOURCE)old/reordr.f \
	$(MESA_HOME)/$(SOURCE)old/trn1e.f \
	$(MESA_HOME)/$(SOURCE)old/trn1so.f \
	$(MESA_HOME)/$(SOURCE)old/trn2e.f \
	$(MESA_HOME)/$(SOURCE)old/totrso.f  \
	$(MESA_HOME)/$(SOURCE)old/vecout.f  


all: $(BINDIR)/$(SOURCE)

$(BINDIR)/$(SOURCE): $(SRCS:.f=.o) $(MESALIB) $(MDLIB)
	$(LOAD_f) -o $(MESA_BIN)/$(SOURCE)

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/$(SOURCE).f
	cat $(SRCS) > ../source/$(SOURCE).f

print: link
	lpr ../source/$(SOURCE).f

clean:
	rm -f *.o core



