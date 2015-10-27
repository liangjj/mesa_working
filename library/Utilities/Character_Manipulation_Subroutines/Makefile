#   @(#)Makefile	5.1 11/6/94
#
#    Makefile for the character library of MesaLib
#    
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR = $(UTILITY)
include $(TOP_DIR)/Makefile.inc
GET = sccs get
REL = 

SRCS = \
	balpar.f \
	chrkey.f \
	crjust.f \
	cskipb.f \
	cskipf.f \
	ctofp.f \
	ctoi.f \
	dcaptl.f \
	ffnext.f \
	fparr.f \
	fpkey.f \
	fptoc.f \
	getfld.f \
	intarr.f \
	intkey.f \
	iskey.f \
	itoc.f \
	keyloc.f \
	keypos.f \
	keystr.f \
	keyval.f \
	lenth.f \
	length.f \
	locase.f \
	logkey.f \
	nxtkey.f \
	nxtlin.f \
	pakstr.f \
	posinp.f \
	positn.f \
	poskey.f \
	putfp.f \
	rmvnb.f \
	skipln.f \
	streqc.f \
	strtyp.f 
OBJECTS = \
	$(SRCS:.f=.o) 

all: $(CHRLIB) 
$(CHRLIB): $(OBJECTS) 
	$(AR) $(ARFLAGS) $(CHRLIB) $(OBJECTS)
	$(RANLIB) $(CHRLIB)
	$(AR) $(ARFLAGS) $(UTILITY_LIB) $(OBJECTS)
	$(RANLIB) $(UTILITY_LIB)
sources: $(SRCS) 
$(SRCS):
	$(GET) $(RELEASE) $@
clean:
	rm -f *.o
