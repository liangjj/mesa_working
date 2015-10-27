#   @(#)Makefile	5.1 11/6/94
#
#    Makefile for the character library of MesaLib
#    
UTILLIB = ../utillib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = fort
FFLAGS = -c -O4
GET = sccs get
REL = 

SRCS = \
	balpar.f \
	chrkey.f \
	chrlen.f \
	crjust.f \
	cskipb.f \
	cskipf.f \
	ctofp.f \
	ctoi.f \
	dgemv.f \
	dger.f \
	dcaptl.f \
	dsyev.f \
	dsbev.f \
	dsterf.f \
	dsteqr.f \
	dlamch.f \
	dlae2.f \
	dlaev2.f \
	dlarf.f \
	dtrmv.f \
	dlansb.f \
	dlansy.f \
	dlapy2.f \
	dlasr.f \
	dlazro.f \
	dstev.f \
	dorgtr.f \
	dsbtrd.f \
	dsytrd.f \
	dlar2v.f \
	dlargv.f \
	dlartg.f \
	dlartv.f \
	dlassq.f \
	dlatrd.f \
	dorgql.f \
	dorgqr.f \
	dsyr2k.f \
	dsytd2.f \
	dlarfb.f \
	dlarfg.f \
	dlarft.f \
	dorgl2.f \
	dnrm2.f \
	dorg2r.f \
	dsymv.f \
	dsyr2.f \
	dtrmm.f \
	drot.f \
	dorg2l.f \
	ffnext.f \
	fparr.f \
	fpkey.f \
	fptoc.f \
	gamfun.f \
	getfld.f \
	intarr.f \
	intkey.f \
	iskey.f \
	itoc.f \
	int_2.f \
	ilaenv.f \
	keyloc.f \
	keypos.f \
	keystr.f \
	keyval.f \
	lsame.f \
	lnkerr.f \
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
	strtyp.f \
	class.f \
	dollar.f \
	gaussq.f \
	gbslve.f \
	gbtql2.f \
	imtql2.f \
	mmove.f \
	prntfm.f \
	prntrm.f \
	pythag.f \
	setlog.f \
	xerbla.f
.f.o:
	$(FC) $(FFLAGS) $<

all: $(UTILLIB)
$(UTILLIB): $(SRCS:.f=.o) 
	$(AR) $(ARFLAGS) $(UTILLIB) $(?:.f=.o)
	$(RANLIB) $(UTILLIB)

sources: $(SRCS) 
$(SRCS):
	$(GET) $(RELEASE) $@

link: sources
	rm -f ../../source/chr.f
	cat $(SRCS) > ../../source/chr.f

print: link
	lpr ../../source/chr.f

clean:
	rm -f *.o
	rm -f core
