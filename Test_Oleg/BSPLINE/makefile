#   make LIBRARY:
 
N  = libBS_08
 
FC     = f90
 
BSRLIBS = ${HOME}/BSR/T/LIBS
BSRMODS = ${HOME}/BSR/T/LIBS

OPTS  = -c -O -static

O  =   spline_param.o           spline_atomic.o          spline_galerkin.o         \
	spline_grid.o            spline_hl.o              spline_integrals.o        \
	spline_moments.o         spline_orbitals.o        spline_slater.o           \
	azl.o                    bav.o                    bcore.o                   \
	bhl.o                    bhwf.o                   bvalue2.o                 \
	bvmv.o                   bxv.o                    bzk.o                     \
	convol.o                 coulom.o                 define_grid.o             \
	define_spline.o          density.o                dinty.o                   \
	dvx.o                    facdyk.o                 facdzk.o                  \
	facsb.o                  facsbl.o                 gauss.o                   \
	grad.o                   hlc.o                    hlm.o                     \
	int_de.o                 int_v.o                  minty.o                   \
	mk.o                     mkc.o                    mkgrid2.o                 \
	mky.o                    mk_moments.o             mmk_cell.o                \
	mmk_diff.o               mnk_cell.o               mnk_diff.o                \
	moments.o                mqk_cell.o               mrk_cell.o                \
	mrk_diff.o               mrm.o                    mtk_cell.o                \
	mtk_diff.o               mvc.o                    mvcv.o                    \
	mvk_cell.o               mvk_diff.o               mwk_cell.o                \
	nk.o                     nkc.o                    nky.o                     \
	nk_moments.o             qk.o                     qkc.o                     \
	qky.o                    qk_moments.o             quadr.o                   \
	rk.o                     rkc.o                    rky.o                     \
	rk_moments.o             r_bwfn.o                 splin3.o                  \
	sum_amb.o                tk.o                     tkc.o                     \
	tky.o                    tk_moments.o             vbsplvd.o                 \
	vinty.o                  vk.o                     vkc.o                     \
	vky.o                    vk_moments.o             wk.o                      \
	wky.o                    wk_moments.o             ykf.o                     \
	yval.o                   zeta_y.o


$(N).a : $(O)
	ar rv $(N).a $(O)
	mv $(N).a $(BSRLIBS)
	mv *.mod $(BSRMODS)
	rm *.o
 
spline_param.o: spline_param.f90
	$(FC) $(OPTS) spline_param.f90
 
spline_atomic.o: spline_atomic.f90
	$(FC) $(OPTS) spline_atomic.f90
 
spline_galerkin.o: spline_galerkin.f90
	$(FC) $(OPTS) spline_galerkin.f90
 
spline_grid.o: spline_grid.f90
	$(FC) $(OPTS) spline_grid.f90
 
spline_hl.o: spline_hl.f90
	$(FC) $(OPTS) spline_hl.f90
 
spline_integrals.o: spline_integrals.f90
	$(FC) $(OPTS) spline_integrals.f90
 
spline_moments.o: spline_moments.f90
	$(FC) $(OPTS) spline_moments.f90
 
spline_orbitals.o: spline_orbitals.f90
	$(FC) $(OPTS) spline_orbitals.f90
 
spline_slater.o: spline_slater.f90
	$(FC) $(OPTS) spline_slater.f90
 
azl.o: azl.f90
	$(FC) $(OPTS) azl.f90
 
bav.o: bav.f90
	$(FC) $(OPTS) bav.f90
 
bcore.o: bcore.f90
	$(FC) $(OPTS) bcore.f90
 
bhl.o: bhl.f90
	$(FC) $(OPTS) bhl.f90
 
bhwf.o: bhwf.f90
	$(FC) $(OPTS) bhwf.f90
 
bvalue2.o: bvalue2.f90
	$(FC) $(OPTS) bvalue2.f90
 
bvmv.o: bvmv.f90
	$(FC) $(OPTS) bvmv.f90
 
bxv.o: bxv.f90
	$(FC) $(OPTS) bxv.f90
 
bzk.o: bzk.f90
	$(FC) $(OPTS) bzk.f90
 
convol.o: convol.f90
	$(FC) $(OPTS) convol.f90
 
coulom.o: coulom.f90
	$(FC) $(OPTS) coulom.f90
 
define_grid.o: define_grid.f90
	$(FC) $(OPTS) define_grid.f90
 
define_spline.o: define_spline.f90
	$(FC) $(OPTS) define_spline.f90
 
density.o: density.f90
	$(FC) $(OPTS) density.f90
 
dinty.o: dinty.f90
	$(FC) $(OPTS) dinty.f90
 
dvx.o: dvx.f90
	$(FC) $(OPTS) dvx.f90
 
facdyk.o: facdyk.f90
	$(FC) $(OPTS) facdyk.f90
 
facdzk.o: facdzk.f90
	$(FC) $(OPTS) facdzk.f90
 
facsb.o: facsb.f90
	$(FC) $(OPTS) facsb.f90
 
facsbl.o: facsbl.f90
	$(FC) $(OPTS) facsbl.f90
 
gauss.o: gauss.f90
	$(FC) $(OPTS) gauss.f90
 
grad.o: grad.f90
	$(FC) $(OPTS) grad.f90
 
hlc.o: hlc.f90
	$(FC) $(OPTS) hlc.f90
 
hlm.o: hlm.f90
	$(FC) $(OPTS) hlm.f90
 
int_de.o: int_de.f90
	$(FC) $(OPTS) int_de.f90
 
int_v.o: int_v.f90
	$(FC) $(OPTS) int_v.f90
 
minty.o: minty.f90
	$(FC) $(OPTS) minty.f90
 
mk.o: mk.f90
	$(FC) $(OPTS) mk.f90
 
mkc.o: mkc.f90
	$(FC) $(OPTS) mkc.f90
 
mkgrid2.o: mkgrid2.f90
	$(FC) $(OPTS) mkgrid2.f90
 
mky.o: mky.f90
	$(FC) $(OPTS) mky.f90
 
mk_moments.o: mk_moments.f90
	$(FC) $(OPTS) mk_moments.f90
 
mmk_cell.o: mmk_cell.f90
	$(FC) $(OPTS) mmk_cell.f90
 
mmk_diff.o: mmk_diff.f90
	$(FC) $(OPTS) mmk_diff.f90
 
mnk_cell.o: mnk_cell.f90
	$(FC) $(OPTS) mnk_cell.f90
 
mnk_diff.o: mnk_diff.f90
	$(FC) $(OPTS) mnk_diff.f90
 
moments.o: moments.f90
	$(FC) $(OPTS) moments.f90
 
mqk_cell.o: mqk_cell.f90
	$(FC) $(OPTS) mqk_cell.f90
 
mrk_cell.o: mrk_cell.f90
	$(FC) $(OPTS) mrk_cell.f90
 
mrk_diff.o: mrk_diff.f90
	$(FC) $(OPTS) mrk_diff.f90
 
mrm.o: mrm.f90
	$(FC) $(OPTS) mrm.f90
 
mtk_cell.o: mtk_cell.f90
	$(FC) $(OPTS) mtk_cell.f90
 
mtk_diff.o: mtk_diff.f90
	$(FC) $(OPTS) mtk_diff.f90
 
mvc.o: mvc.f90
	$(FC) $(OPTS) mvc.f90
 
mvcv.o: mvcv.f90
	$(FC) $(OPTS) mvcv.f90
 
mvk_cell.o: mvk_cell.f90
	$(FC) $(OPTS) mvk_cell.f90
 
mvk_diff.o: mvk_diff.f90
	$(FC) $(OPTS) mvk_diff.f90
 
mwk_cell.o: mwk_cell.f90
	$(FC) $(OPTS) mwk_cell.f90
 
nk.o: nk.f90
	$(FC) $(OPTS) nk.f90
 
nkc.o: nkc.f90
	$(FC) $(OPTS) nkc.f90
 
nky.o: nky.f90
	$(FC) $(OPTS) nky.f90
 
nk_moments.o: nk_moments.f90
	$(FC) $(OPTS) nk_moments.f90
 
qk.o: qk.f90
	$(FC) $(OPTS) qk.f90
 
qkc.o: qkc.f90
	$(FC) $(OPTS) qkc.f90
 
qky.o: qky.f90
	$(FC) $(OPTS) qky.f90
 
qk_moments.o: qk_moments.f90
	$(FC) $(OPTS) qk_moments.f90
 
quadr.o: quadr.f90
	$(FC) $(OPTS) quadr.f90
 
rk.o: rk.f90
	$(FC) $(OPTS) rk.f90
 
rkc.o: rkc.f90
	$(FC) $(OPTS) rkc.f90
 
rky.o: rky.f90
	$(FC) $(OPTS) rky.f90
 
rk_moments.o: rk_moments.f90
	$(FC) $(OPTS) rk_moments.f90
 
r_bwfn.o: r_bwfn.f90
	$(FC) $(OPTS) r_bwfn.f90
 
splin3.o: splin3.f90
	$(FC) $(OPTS) splin3.f90
 
sum_amb.o: sum_amb.f90
	$(FC) $(OPTS) sum_amb.f90
 
tk.o: tk.f90
	$(FC) $(OPTS) tk.f90
 
tkc.o: tkc.f90
	$(FC) $(OPTS) tkc.f90
 
tky.o: tky.f90
	$(FC) $(OPTS) tky.f90
 
tk_moments.o: tk_moments.f90
	$(FC) $(OPTS) tk_moments.f90
 
vbsplvd.o: vbsplvd.f90
	$(FC) $(OPTS) vbsplvd.f90
 
vinty.o: vinty.f90
	$(FC) $(OPTS) vinty.f90
 
vk.o: vk.f90
	$(FC) $(OPTS) vk.f90
 
vkc.o: vkc.f90
	$(FC) $(OPTS) vkc.f90
 
vky.o: vky.f90
	$(FC) $(OPTS) vky.f90
 
vk_moments.o: vk_moments.f90
	$(FC) $(OPTS) vk_moments.f90
 
wk.o: wk.f90
	$(FC) $(OPTS) wk.f90
 
wky.o: wky.f90
	$(FC) $(OPTS) wky.f90
 
wk_moments.o: wk_moments.f90
	$(FC) $(OPTS) wk_moments.f90
 
ykf.o: ykf.f90
	$(FC) $(OPTS) ykf.f90
 
yval.o: yval.f90
	$(FC) $(OPTS) yval.f90
 
zeta_y.o: zeta_y.f90
	$(FC) $(OPTS) zeta_y.f90
 
 
