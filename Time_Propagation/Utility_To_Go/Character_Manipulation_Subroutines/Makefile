#   @(#)Makefile	5.1 11/6/94
#
#    Makefile for the character library of Util_Lib
#    
RANLIB = ranlib

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
	length.f \
	lenth.f \
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
 
.f.o:
	$(F90) $(F90FLAGS) $<

all: $(UTILITY_LIB)
$(UTILITY_LIB): $(SRCS:.f=.o) 
	$(AR) $(ARFLAGS) $(UTILITY_LIB) $(?:.f=.o)
	$(RANLIB) $(UTILITY_LIB)

sources: $(SRCS) 
$(SRCS):
	$(GET) $(RELEASE) $@

clean:
	rm -f *.o
