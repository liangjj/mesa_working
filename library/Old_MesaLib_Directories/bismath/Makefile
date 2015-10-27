#    %W%   %G%
#
#    Makefile for the utility portion of Mylib
#    
MYLIB = ../mylib.a
MESALIB = ../mesalib.a
AR = -ar
RANLIB = ranlib
FFLAGS = -c
REL = 

SRCS = \
	asub.f \
	banmat.f \
	bessel.f \
	besj0.f \
        besj1.f \
        besy0.f \
        besy1.f \
	besjy01.f \
	besjyn.f \
	binom.f \
	bffgh.f \
	bser.f \
	bsplin.f \
	bsplpp.f \
	bsplvb.f \
	cot.f \
	cgvmmul.f \
	chnvar.f \
	class.f \
	classp.f \
	cmadd.f \
	cpoly.f \
	csmul.f \
	cvfill.f \
	cvmmul.f \
	cvmul.f \
	cvscal.f \
	d9b0mp.f \
	d9b1mp.f \
	envj.f \
	fact.f \
	fndbrk.f \
	fzero.f \
	gamfun.f \
	gasdev.f \
	gaussp.f \
	gaussq.f \
	gbslve.f \
	gbtql2.f \
	getqpt.f \
	gpoly.f \
	ifill.f \
	init.f \
	initds.f \
	int_2.f \
	int_3.f \
	int_4.f \
	int_5.f \
	int_6.f \
	int_7.f \
	int_8.f \
	int_9.f \
	lancz.f \
	legend.f \
	lgngr2.f \
	lgngr.f \
	lgpoly.f \
	lgndx2.f \
	lgrply.f \
	lpoly.f \
	madd.f \
	mmove.f \
	mscale.f \
	mvmul.f \
	mvscal.f \
	msta1.f \
	msta2.f \
	modbes.f \
	nwtrap.f \
	plyfit.f \
	polfit.f \
	polvlu.f \
	poly0.f \
	polyab.f \
	ppval.f \
	ppvalu.f \
	prespl.f \
	ran1.f \
	rbes.f \
        rcbes.f \
        rc1bes.f \
	rcbesb.f \
	rcbesf.f \
	reight.f \
	rfive.f \
	rfour.f \
	rseven.f \
	rsix.f \
	rthree.f \
	rtwo.f \
	rnine.f \
	secder.f \
	sgefs.f \
	smuls.f \
	snorm.f \
	splcof.f \
	splmat.f \
	sumel.f \
	udagv.f \
	udagmu.f \
	v2c.f \
	v2r.f \
	vabs.f \
	vec2di.f \
	vimmul.f \
	vmmul.f \
	vscale.f \
	vwxys.f


.f.o:
	$(FC) $(FFLAGS) $<

all: $(MYLIB)
$(MYLIB): $(SRCS:.f=.o)
	$(AR) rvu $(MYLIB) $(?:.f=.o)
	$(RANLIB) $(MYLIB)

all: $(MESALIB)
$(MESALIB): $(SRCS:.f=.o)
	$(AR) rvu $(MESALIB) $(?:.f=.o)
	$(RANLIB) $(MESALIB)
