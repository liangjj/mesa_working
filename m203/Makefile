#   @(#)Makefile	5.1   11/6/94
#
# Makefile for M101
#
include $(MESA_HOME)/Makefile.inc
SOURCE=m203
GET = sccs get
REL = 

SRCS = \
	$(MESA_HOME)/$(SOURCE)/cztran.f \
	$(MESA_HOME)/$(SOURCE)/eigint.f \
	$(MESA_HOME)/$(SOURCE)/fillop.f \
	$(MESA_HOME)/$(SOURCE)/$(SOURCE).f \
	$(MESA_HOME)/$(SOURCE)/schmidt.f \
	$(MESA_HOME)/$(SOURCE)/trmain.f \
	$(MESA_HOME)/$(SOURCE)/tvec.f

all: $(MESA_BIN)/$(SOURCE) 

$(MESA_BIN)/$(SOURCE): $(SRCS:.f=.o)
	$(LOAD) -o $(MESA_BIN)/$(SOURCE)

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f $(MESA_HOME)/source/$(SOURCE).f
	cat $(SRCS) > $(MESA_HOME)/source/$(SOURCE).f

print: link
	lpr $(MESA_HOME)/source/$(SOURCE).f

clean:
	rm -f *.o core
