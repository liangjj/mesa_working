#    @(#)Makefile	5.2   4/17/95
#
#    Makefile for the mathematics library of MesaLib
#    
.SUFFIXES: .f90 .mod .F90 $(SUFFIXES)
TOP_DIR = $(UTILITY)
include $(TOP_DIR)/Makefile.inc
#
SRCS_MESA = \
	ambc.f \
	ambct.f \
	ambtc.f \
	ambtct.f \
	apbc.f \
	apbct.f \
	apbctt.f \
	apbtc.f \
	apbtct.f \
	arrmax.f \
	atebc.f \
	atembc.f \
	atmbc.f \
	atpbc.f \
	atpbct.f \
	bico.f \
	bliu.f \
	coulomb.f \
	degrsp.f \
	ebc.f \
	ebct.f \
	ebtc.f \
	ebtct.f \
	embc.f \
	embct.f \
	embtc.f \
	embtct.f \
	factln.f \
	factrl.f \
	flin.f \
	fold.f \
	gammln.f \
	givens.f \
	gvadd.f \
	gvdiv.f \
	gvmul.f \
	gvsub.f \
	iveceq.f \
	izero.f \
	lowdin.f \
	minvrt.f \
	runit.f \
	rveceq.f \
	rzero.f \
	sadd.f \
	schmdt.f \
	schmidt.f \
	sgmm.f \
	sinv.f \
	smul.f \
	solver.f \
	sorter.f \
	sorter1.f \
	sorter2.f \
	sorter3.f \
	sorter4.f \
	sqtotr.f \
	ssub.f \
	trace.f \
	trtosq.f \
	vadd.f \
	vclean.f \
	vdiv.f \
	vexp.f \
	vfill.f \
	vinv.f \
	vmove.f \
	vmul.f \
	vadmul.f \
	vneg.f \
	vpoly.f \
	vpolyd.f \
	vpower.f \
	vsqrt.f \
	vsub.f \
	vwx.f \
	vwxs.f \
	vwxy.f \
	vwxyi.f 

SRCS_BIS = \
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

SRCS_BISf90 = \
	F_3J.f90 \
	Z_3J.f90 \
	Z_6J.f90 \
	Z_9J.f90 \
	Z_6JJ.f90 

OBJECT.f = \
	$(SRCS_MESA:.f=.o) \
	$(SRCS_BIS:.f=.o) 

OBJECT.f90 = \
	$(SRCS_BISf90:.f90=.o) 

all: $(MATHLIB) $(OBJECT.f) $(OBJECT.f90)
$(MATHLIB): $(OBJECT.f) $(OBJECT.f90)
	$(AR) $(ARFLAGS) $(MATHLIB) $(OBJECT.f)
	$(RANLIB) $(MATHLIB)
	$(AR) $(ARFLAGS) $(MATHLIB) $(OBJECT.f90)
	$(RANLIB) $(MATHLIB)
	$(AR) $(ARFLAGS) $(UTILITY_LIB) $(OBJECT.f)
	$(RANLIB) $(UTILITY_LIB)
	$(AR) $(ARFLAGS) $(UTILITY_LIB) $(OBJECT.f90)
	$(RANLIB) $(UTILITY_LIB)

clean:
	rm -f *.o

