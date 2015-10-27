#   @(#)Makefile        1.4   8/2/91
#
# Here too a $(BSPLINE_SOURCE_LIB) is produced which needs to be loades with any code using the BSPLINE program
#
TOP_DIR = $(MESA_HOME)/Test_Oleg
include $(TOP_DIR)/Makefile.inc
SRCS = \
	azl.f90 \
	bav.f90 \
	bcore.f90 \
	bhl.f90 \
	bhwf.f90 \
	bvalue2.f90 \
	bvmv.f90 \
	bxv.f90 \
	bzk.f90 \
	convol.f90 \
	coulom.f90 \
	define_grid.f90 \
	define_spline.f90 \
	density.f90 \
	dinty.f90 \
	dvx.f90 \
	facdyk.f90 \
	facdzk.f90 \
	facsb.f90 \
	facsbl.f90 \
	gauss.f90 \
	grad.f90 \
	hlc.f90 \
	hlm.f90 \
	int_de.f90 \
	int_v.f90 \
	minty.f90 \
	mk.f90 \
	mkc.f90 \
	mkgrid2.f90 \
	mky.f90 \
	mk_moments.f90 \
	mmk_cell.f90 \
	mmk_diff.f90 \
	mnk_cell.f90 \
	mnk_diff.f90 \
	moments.f90 \
	mqk_cell.f90 \
	mrk_cell.f90 \
	mrk_diff.f90 \
	mrm.f90 \
	mtk_cell.f90 \
	mtk_diff.f90 \
	mvc.f90 \
	mvcv.f90 \
	mvk_cell.f90 \
	mvk_diff.f90 \
	mwk_cell.f90 \
	nk.f90 \
	nkc.f90 \
	nky.f90 \
	nk_moments.f90 \
	qk.f90 \
	qkc.f90 \
	qky.f90 \
	qk_moments.f90 \
	quadr.f90 \
	rk.f90 \
	rkc.f90 \
	rky.f90 \
	rk_moments.f90 \
	r_bwfn.f90 \
	splin3.f90 \
	sum_amb.f90 \
	tk.f90 \
	tkc.f90 \
	tky.f90 \
	tk_moments.f90 \
	vbsplvd.f90 \
	vinty.f90 \
	vk.f90 \
	vkc.f90 \
	vky.f90 \
	vk_moments.f90 \
	wk.f90 \
	wky.f90 \
	wk_moments.f90 \
	ykf.f90 \
	yval.f90 \
	zeta_y.f90

all: $(BSPLINE_SOURCE_LIB)
$(BSPLINE_SOURCE_LIB): $(SRCS:.f90=.o)
	$(AR) $(ARFLAGS) $(BSPLINE_SOURCE_LIB) $(?:.f90=.o)
	$(RANLIB) $(BSPLINE_SOURCE_LIB)
clean:
	rm -f *.o
	rm -f *.mod

