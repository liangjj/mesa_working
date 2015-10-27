#    @(#)Makefile	5.1   11/6/94
#
#    Makefile for the utility portion of MesaLib
#    
UTIL_LIB = ../util_lib_f90.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
GET = sccs get
REL = 

SRCS_MESA = \
	basis.f \
	chainx.f \
	chkpnt.f \
	corpr1.f \
	dismat.f \
	dollar.f \
	drum.f \
	enforc.f \
	filchg.f \
	fillel.f \
	filmas.f \
	fixdrt.f \
	fmstrt.f \
	frqprt.f \
	funcnm.f \
	getcf.f \
	liopen.f \
	lnkerr.f \
	masswt.f \
	matout.f \
	matprt.f \
	norm1.f \
	normal.f \
	print.f \
	project.f \
	prsrte.f \
	rcoord.f \
	replic.f \
	rtelbl.f \
	rzmat.f \
	rzsub.f \
	subvar.f \
	togenc.f \
	vec.f \
	vecout.f \
	vprod.f \
	wcoord.f \
	wlmat.f \
	wmat.f \
	wvec.f \
	wzmat.f \
	wzsub.f \
	zconv.f \
	zprint.f \
	ztoc.f

SRCS_BIS = \
	abschm.f \
	aeqbc.f \
	ambcc.f \
	ambcct.f \
	ambccx.f \
	ambcxx.f \
	ambtcc.f \
	ambtccx.f \
	ambtcx.f \
	ambtcxx.f \
	ambctxx.f \
	amcbc.f \
	amcbct.f \
	amcbctx.f \
	amcbcx.f \
	amcbtc.f \
	amcbtcx.f \
	apbcc.f \
	apbccx.f \
	apbcxx.f \
	apbtcxx.f \
	apbctxx.f \
	apbtcc.f \
	apbtccx.f \
	apcbc.f \
	apcbct.f \
	apcbctx.f \
	apcbtc.f \
	apcbtcx.f \
	bfunct.f \
	cambc.f \
	cambct.f \
	cambctx.f \
	cambcx.f \
	cambtc.f \
	cambtcx.f \
	capbc.f \
	capbct.f \
	capbctx.f \
	capbcx.f \
	capbtc.f \
	capbtcx.f \
	cardin.f \
	cc2opy.f \
	cebc.f \
	cebct.f \
	cebctx.f \
	cebcx.f \
	cebtc.f \
	cehbtc.f \
	cebhct.f \
	cehbtcx.f \
	cebtct.f \
	cebtcx.f \
	cmprir.f \
	copy.f \
	crvadd.f \
	cschmt.f \
	csxpy.f \
	cvadd.f \
	cvconj.f \
	czero.f \
	dfactl.f \
	dagfil.f \
	dschmt.f \
	ebcc.f \
	ebccx.f \
	ebctx.f \
	ebctxx.f \
	ebcx.f \
	ebcxx.f \
	ebtcc.f \
	ebtccx.f \
	ebtcx.f \
	ebtcxx.f \
	ecbc.f \
	ecbct.f \
	ecbcx.f \
	ecbtc.f \
	erf.f \
	erfc.f \
	etimes.f \
	gam0.f \
	gfunct.f \
	gramsc.f \
	gschmt.f \
	factl.f \
        filtri.f \
	hamslv.f \
	icopy.f \
	lgndrx.f \
	linslv.f \
	mprii.f \
	mprir.f \
	mprri.f \
	mprrr.f \
	order.f \
	prntcm.f \
	prntcmn.f \
	prntfmn.f \
	prntfm.f \
	prntrm.f \
	rtocm.f \
	schmab.f \
	series.f \
	setlog.f \
	setprn.f \
	trimm.f \
	trimv.f \
	va2apbc.f \
	wrcvec.f \
        wrvec.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(UTIL_LIB) 
$(UTIL_LIB): $(SRCS_MESA:.f=.o) $(SRCS_BIS:.f=.o)
	$(AR) $(ARFLAGS) $(UTIL_LIB) $(?:.f=.o)
	$(RANLIB) $(UTIL_LIB)

sources: $(SRCS_MESA) $(SRCS_BIS)
$(SRCS_MESA) $(SRCS_BIS):
	$(GET) $(RELEASE) $@

clean:
	rm -f *.o
