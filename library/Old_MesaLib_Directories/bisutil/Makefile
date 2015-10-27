#   
#
#    Makefile for the utility portion of Mesctlib
#    
MYLIB = ../mylib.a
MESALIB = ../mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib
FFLAGS = -c
REL = 

SRCS = \
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
	prntcmn.f \
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
	$(FC) $(FFLAGS)  $<

all: $(MYLIB) 
$(MYLIB): $(SRCS:.f=.o)
	$(AR) rvu $(MYLIB) $(?:.f=.o)
	$(RANLIB) $(MYLIB)

all: $(MESALIB) 
$(MESALIB): $(SRCS:.f=.o)
	$(AR) rvu $(MESALIB) $(?:.f=.o)
	$(RANLIB) $(MESALIB)
