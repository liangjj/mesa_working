# Microsoft Developer Studio Generated NMAKE File, Based on BS.dsp
!IF "$(CFG)" == ""
CFG=BS - WIN32 RELEASE
!MESSAGE No configuration specified. Defaulting to BS - WIN32 RELEASE.
!ENDIF 

!IF "$(CFG)" != "BS - Win32 Release"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "BS.mak" CFG="BS - WIN32 RELEASE"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "BS - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

CPP=cl.exe
F90=df.exe
RSC=rc.exe
OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "..\..\..\LIB\BS.lib" "$(OUTDIR)\spline_atomic.mod" "$(OUTDIR)\spline_galerkin.mod" "$(OUTDIR)\spline_grid.mod" "$(OUTDIR)\spline_hl.mod" "$(OUTDIR)\spline_integrals.mod" "$(OUTDIR)\spline_moments.mod" "$(OUTDIR)\spline_orbitals.mod" "$(OUTDIR)\spline_param.mod" "$(OUTDIR)\spline_slater.mod" "$(OUTDIR)\spline_Rk_integrals.mod" "$(OUTDIR)\spline_densities.mod"


CLEAN :
	-@erase "$(INTDIR)\azl.obj"
	-@erase "$(INTDIR)\bav.obj"
	-@erase "$(INTDIR)\bcore.obj"
	-@erase "$(INTDIR)\bhl.obj"
	-@erase "$(INTDIR)\bhwf.obj"
	-@erase "$(INTDIR)\bvalue2.obj"
	-@erase "$(INTDIR)\bvmv.obj"
	-@erase "$(INTDIR)\bxv.obj"
	-@erase "$(INTDIR)\bzk.obj"
	-@erase "$(INTDIR)\convol.obj"
	-@erase "$(INTDIR)\coulom.obj"
	-@erase "$(INTDIR)\define_grid.obj"
	-@erase "$(INTDIR)\define_spline.obj"
	-@erase "$(INTDIR)\density.obj"
	-@erase "$(INTDIR)\dinty.obj"
	-@erase "$(INTDIR)\dvx.obj"
	-@erase "$(INTDIR)\facdyk.obj"
	-@erase "$(INTDIR)\facdzk.obj"
	-@erase "$(INTDIR)\facsb.obj"
	-@erase "$(INTDIR)\facsbl.obj"
	-@erase "$(INTDIR)\gauss.obj"
	-@erase "$(INTDIR)\grad.obj"
	-@erase "$(INTDIR)\hlc.obj"
	-@erase "$(INTDIR)\hlm.obj"
	-@erase "$(INTDIR)\int_de.obj"
	-@erase "$(INTDIR)\int_v.obj"
	-@erase "$(INTDIR)\minty.obj"
	-@erase "$(INTDIR)\mk.obj"
	-@erase "$(INTDIR)\mk_moments.obj"
	-@erase "$(INTDIR)\mkc.obj"
	-@erase "$(INTDIR)\mky.obj"
	-@erase "$(INTDIR)\mmk_cell.obj"
	-@erase "$(INTDIR)\mmk_diff.obj"
	-@erase "$(INTDIR)\mnk_cell.obj"
	-@erase "$(INTDIR)\mnk_diff.obj"
	-@erase "$(INTDIR)\moments.obj"
	-@erase "$(INTDIR)\mqk_cell.obj"
	-@erase "$(INTDIR)\mrk_cell.obj"
	-@erase "$(INTDIR)\mrk_diff.obj"
	-@erase "$(INTDIR)\mrm.obj"
	-@erase "$(INTDIR)\mtk_cell.obj"
	-@erase "$(INTDIR)\mtk_diff.obj"
	-@erase "$(INTDIR)\mvc.obj"
	-@erase "$(INTDIR)\mvcv.obj"
	-@erase "$(INTDIR)\mvk_cell.obj"
	-@erase "$(INTDIR)\mvk_diff.obj"
	-@erase "$(INTDIR)\mwk_cell.obj"
	-@erase "$(INTDIR)\nk.obj"
	-@erase "$(INTDIR)\nk_moments.obj"
	-@erase "$(INTDIR)\nkc.obj"
	-@erase "$(INTDIR)\nky.obj"
	-@erase "$(INTDIR)\qk.obj"
	-@erase "$(INTDIR)\qk_moments.obj"
	-@erase "$(INTDIR)\qkc.obj"
	-@erase "$(INTDIR)\qky.obj"
	-@erase "$(INTDIR)\quadr.obj"
	-@erase "$(INTDIR)\r_bwfn.obj"
	-@erase "$(INTDIR)\rk.obj"
	-@erase "$(INTDIR)\rk_moments.obj"
	-@erase "$(INTDIR)\rkc.obj"
	-@erase "$(INTDIR)\rky.obj"
	-@erase "$(INTDIR)\splin3.obj"
	-@erase "$(INTDIR)\spline_atomic.mod"
	-@erase "$(INTDIR)\spline_atomic.obj"
	-@erase "$(INTDIR)\spline_densities.mod"
	-@erase "$(INTDIR)\spline_densities.obj"
	-@erase "$(INTDIR)\spline_galerkin.mod"
	-@erase "$(INTDIR)\spline_galerkin.obj"
	-@erase "$(INTDIR)\spline_grid.mod"
	-@erase "$(INTDIR)\spline_grid.obj"
	-@erase "$(INTDIR)\spline_hl.mod"
	-@erase "$(INTDIR)\spline_hl.obj"
	-@erase "$(INTDIR)\spline_integrals.mod"
	-@erase "$(INTDIR)\spline_integrals.obj"
	-@erase "$(INTDIR)\spline_moments.mod"
	-@erase "$(INTDIR)\spline_moments.obj"
	-@erase "$(INTDIR)\spline_orbitals.mod"
	-@erase "$(INTDIR)\spline_orbitals.obj"
	-@erase "$(INTDIR)\spline_param.mod"
	-@erase "$(INTDIR)\spline_param.obj"
	-@erase "$(INTDIR)\spline_Rk_integrals.mod"
	-@erase "$(INTDIR)\spline_Rk_integrals.obj"
	-@erase "$(INTDIR)\spline_slater.mod"
	-@erase "$(INTDIR)\spline_slater.obj"
	-@erase "$(INTDIR)\sum_amb.obj"
	-@erase "$(INTDIR)\tk.obj"
	-@erase "$(INTDIR)\tk_moments.obj"
	-@erase "$(INTDIR)\tkc.obj"
	-@erase "$(INTDIR)\tky.obj"
	-@erase "$(INTDIR)\vbsplvd.obj"
	-@erase "$(INTDIR)\vinty.obj"
	-@erase "$(INTDIR)\vk.obj"
	-@erase "$(INTDIR)\vk_moments.obj"
	-@erase "$(INTDIR)\vkc.obj"
	-@erase "$(INTDIR)\vky.obj"
	-@erase "$(INTDIR)\wk.obj"
	-@erase "$(INTDIR)\wk_moments.obj"
	-@erase "$(INTDIR)\wky.obj"
	-@erase "$(INTDIR)\ykf.obj"
	-@erase "$(INTDIR)\yval.obj"
	-@erase "$(INTDIR)\zeta_y.obj"
	-@erase "..\..\..\LIB\BS.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\BS.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"C:\oleg\LIB\BS.lib" 
LIB32_OBJS= \
	"$(INTDIR)\azl.obj" \
	"$(INTDIR)\bav.obj" \
	"$(INTDIR)\bcore.obj" \
	"$(INTDIR)\bhl.obj" \
	"$(INTDIR)\bhwf.obj" \
	"$(INTDIR)\bvalue2.obj" \
	"$(INTDIR)\bvmv.obj" \
	"$(INTDIR)\bxv.obj" \
	"$(INTDIR)\bzk.obj" \
	"$(INTDIR)\convol.obj" \
	"$(INTDIR)\coulom.obj" \
	"$(INTDIR)\define_grid.obj" \
	"$(INTDIR)\define_spline.obj" \
	"$(INTDIR)\density.obj" \
	"$(INTDIR)\dinty.obj" \
	"$(INTDIR)\dvx.obj" \
	"$(INTDIR)\facdyk.obj" \
	"$(INTDIR)\facdzk.obj" \
	"$(INTDIR)\facsb.obj" \
	"$(INTDIR)\facsbl.obj" \
	"$(INTDIR)\gauss.obj" \
	"$(INTDIR)\grad.obj" \
	"$(INTDIR)\hlc.obj" \
	"$(INTDIR)\hlm.obj" \
	"$(INTDIR)\int_de.obj" \
	"$(INTDIR)\int_v.obj" \
	"$(INTDIR)\minty.obj" \
	"$(INTDIR)\mk.obj" \
	"$(INTDIR)\mk_moments.obj" \
	"$(INTDIR)\mkc.obj" \
	"$(INTDIR)\mky.obj" \
	"$(INTDIR)\mmk_cell.obj" \
	"$(INTDIR)\mmk_diff.obj" \
	"$(INTDIR)\mnk_cell.obj" \
	"$(INTDIR)\mnk_diff.obj" \
	"$(INTDIR)\moments.obj" \
	"$(INTDIR)\mqk_cell.obj" \
	"$(INTDIR)\mrk_cell.obj" \
	"$(INTDIR)\mrk_diff.obj" \
	"$(INTDIR)\mrm.obj" \
	"$(INTDIR)\mtk_cell.obj" \
	"$(INTDIR)\mtk_diff.obj" \
	"$(INTDIR)\mvc.obj" \
	"$(INTDIR)\mvcv.obj" \
	"$(INTDIR)\mvk_cell.obj" \
	"$(INTDIR)\mvk_diff.obj" \
	"$(INTDIR)\mwk_cell.obj" \
	"$(INTDIR)\nk.obj" \
	"$(INTDIR)\nk_moments.obj" \
	"$(INTDIR)\nkc.obj" \
	"$(INTDIR)\nky.obj" \
	"$(INTDIR)\qk.obj" \
	"$(INTDIR)\qk_moments.obj" \
	"$(INTDIR)\qkc.obj" \
	"$(INTDIR)\qky.obj" \
	"$(INTDIR)\quadr.obj" \
	"$(INTDIR)\r_bwfn.obj" \
	"$(INTDIR)\rk.obj" \
	"$(INTDIR)\rk_moments.obj" \
	"$(INTDIR)\rkc.obj" \
	"$(INTDIR)\rky.obj" \
	"$(INTDIR)\splin3.obj" \
	"$(INTDIR)\spline_atomic.obj" \
	"$(INTDIR)\spline_galerkin.obj" \
	"$(INTDIR)\spline_grid.obj" \
	"$(INTDIR)\spline_hl.obj" \
	"$(INTDIR)\spline_integrals.obj" \
	"$(INTDIR)\spline_moments.obj" \
	"$(INTDIR)\spline_orbitals.obj" \
	"$(INTDIR)\spline_param.obj" \
	"$(INTDIR)\spline_slater.obj" \
	"$(INTDIR)\sum_amb.obj" \
	"$(INTDIR)\tk.obj" \
	"$(INTDIR)\tk_moments.obj" \
	"$(INTDIR)\tkc.obj" \
	"$(INTDIR)\tky.obj" \
	"$(INTDIR)\vbsplvd.obj" \
	"$(INTDIR)\vinty.obj" \
	"$(INTDIR)\vk.obj" \
	"$(INTDIR)\vk_moments.obj" \
	"$(INTDIR)\vkc.obj" \
	"$(INTDIR)\vky.obj" \
	"$(INTDIR)\wk.obj" \
	"$(INTDIR)\wk_moments.obj" \
	"$(INTDIR)\wky.obj" \
	"$(INTDIR)\ykf.obj" \
	"$(INTDIR)\yval.obj" \
	"$(INTDIR)\zeta_y.obj" \
	"$(INTDIR)\spline_densities.obj" \
	"$(INTDIR)\spline_Rk_integrals.obj"

"..\..\..\LIB\BS.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

SOURCE="$(InputPath)"
DS_POSTBUILD_DEP=$(INTDIR)\postbld.dep

ALL : $(DS_POSTBUILD_DEP)

# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

$(DS_POSTBUILD_DEP) : "..\..\..\LIB\BS.lib" "$(OUTDIR)\spline_atomic.mod" "$(OUTDIR)\spline_galerkin.mod" "$(OUTDIR)\spline_grid.mod" "$(OUTDIR)\spline_hl.mod" "$(OUTDIR)\spline_integrals.mod" "$(OUTDIR)\spline_moments.mod" "$(OUTDIR)\spline_orbitals.mod" "$(OUTDIR)\spline_param.mod" "$(OUTDIR)\spline_slater.mod" "$(OUTDIR)\spline_Rk_integrals.mod" "$(OUTDIR)\spline_densities.mod"
   copy release\*.mod C:\oleg\INCLUDE
	copy release\*.mod .
	echo Helper for Post-build step > "$(DS_POSTBUILD_DEP)"

CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /Fp"$(INTDIR)\BS.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

F90_PROJ=/compile_only /include:"$(INTDIR)\\" /nologo /warn:argument_checking /warn:declarations /warn:nofileopt /module:"Release/" /object:"Release/" 

.SUFFIXES: .fpp

F90_OBJS=.\Release/

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("BS.dep")
!INCLUDE "BS.dep"
!ELSE 
!MESSAGE Warning: cannot find "BS.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "BS - Win32 Release"
SOURCE=.\azl.f90

"$(INTDIR)\azl.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\bav.f90

"$(INTDIR)\bav.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\bcore.f90

"$(INTDIR)\bcore.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\bhl.f90

"$(INTDIR)\bhl.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\bhwf.f90

"$(INTDIR)\bhwf.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\bvalue2.f90

"$(INTDIR)\bvalue2.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\bvmv.f90

"$(INTDIR)\bvmv.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\bxv.f90

"$(INTDIR)\bxv.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\bzk.f90

"$(INTDIR)\bzk.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\convol.f90

"$(INTDIR)\convol.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\coulom.f90

"$(INTDIR)\coulom.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\define_grid.f90

"$(INTDIR)\define_grid.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\define_spline.f90

"$(INTDIR)\define_spline.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\density.f90

"$(INTDIR)\density.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\dinty.f90

"$(INTDIR)\dinty.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\dvx.f90

"$(INTDIR)\dvx.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\facdyk.f90

"$(INTDIR)\facdyk.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\facdzk.f90

"$(INTDIR)\facdzk.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\facsb.f90

"$(INTDIR)\facsb.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\facsbl.f90

"$(INTDIR)\facsbl.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\gauss.f90

"$(INTDIR)\gauss.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\grad.f90

"$(INTDIR)\grad.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\hlc.f90

"$(INTDIR)\hlc.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\hlm.f90

"$(INTDIR)\hlm.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\int_de.f90

"$(INTDIR)\int_de.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\int_v.f90

"$(INTDIR)\int_v.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\minty.f90

"$(INTDIR)\minty.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mk.f90

"$(INTDIR)\mk.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mk_moments.f90

"$(INTDIR)\mk_moments.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mkc.f90

"$(INTDIR)\mkc.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mky.f90

"$(INTDIR)\mky.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mmk_cell.f90

"$(INTDIR)\mmk_cell.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mmk_diff.f90

"$(INTDIR)\mmk_diff.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mnk_cell.f90

"$(INTDIR)\mnk_cell.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mnk_diff.f90

"$(INTDIR)\mnk_diff.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\moments.f90

"$(INTDIR)\moments.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mqk_cell.f90

"$(INTDIR)\mqk_cell.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mrk_cell.f90

"$(INTDIR)\mrk_cell.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mrk_diff.f90

"$(INTDIR)\mrk_diff.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mrm.f90

"$(INTDIR)\mrm.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mtk_cell.f90

"$(INTDIR)\mtk_cell.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mtk_diff.f90

"$(INTDIR)\mtk_diff.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mvc.f90

"$(INTDIR)\mvc.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mvcv.f90

"$(INTDIR)\mvcv.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mvk_cell.f90

"$(INTDIR)\mvk_cell.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mvk_diff.f90

"$(INTDIR)\mvk_diff.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\mwk_cell.f90

"$(INTDIR)\mwk_cell.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\nk.f90

"$(INTDIR)\nk.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\nk_moments.f90

"$(INTDIR)\nk_moments.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\nkc.f90

"$(INTDIR)\nkc.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\nky.f90

"$(INTDIR)\nky.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\qk.f90

"$(INTDIR)\qk.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\qk_moments.f90

"$(INTDIR)\qk_moments.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\qkc.f90

"$(INTDIR)\qkc.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\qky.f90

"$(INTDIR)\qky.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\quadr.f90

"$(INTDIR)\quadr.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\r_bwfn.f90

"$(INTDIR)\r_bwfn.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\rk.f90

"$(INTDIR)\rk.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\rk_moments.f90

"$(INTDIR)\rk_moments.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\rkc.f90

"$(INTDIR)\rkc.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\rky.f90

"$(INTDIR)\rky.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\splin3.f90

"$(INTDIR)\splin3.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\spline_atomic.f90
F90_MODOUT=\
	"spline_atomic"


"$(INTDIR)\spline_atomic.obj"	"$(INTDIR)\spline_atomic.mod" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\spline_densities.f90
F90_MODOUT=\
	"spline_densities"


"$(INTDIR)\spline_densities.obj"	"$(INTDIR)\spline_densities.mod" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\spline_galerkin.f90
F90_MODOUT=\
	"spline_galerkin"


"$(INTDIR)\spline_galerkin.obj"	"$(INTDIR)\spline_galerkin.mod" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\spline_grid.f90
F90_MODOUT=\
	"spline_grid"


"$(INTDIR)\spline_grid.obj"	"$(INTDIR)\spline_grid.mod" : $(SOURCE) "$(INTDIR)"
	$(F90) /assume:nosource_include /compile_only /include:"$(INTDIR)\\" /nologo /warn:argument_checking /warn:declarations /warn:nofileopt /module:"Release/" /object:"Release/" $(SOURCE)


SOURCE=.\spline_hl.f90
F90_MODOUT=\
	"spline_hl"


"$(INTDIR)\spline_hl.obj"	"$(INTDIR)\spline_hl.mod" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\spline_integrals.f90
F90_MODOUT=\
	"spline_integrals"


"$(INTDIR)\spline_integrals.obj"	"$(INTDIR)\spline_integrals.mod" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\spline_moments.f90
F90_MODOUT=\
	"spline_moments"


"$(INTDIR)\spline_moments.obj"	"$(INTDIR)\spline_moments.mod" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\spline_orbitals.f90
F90_MODOUT=\
	"spline_orbitals"


"$(INTDIR)\spline_orbitals.obj"	"$(INTDIR)\spline_orbitals.mod" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\spline_param.f90
F90_MODOUT=\
	"spline_param"


"$(INTDIR)\spline_param.obj"	"$(INTDIR)\spline_param.mod" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\spline_Rk_integrals.f90
F90_MODOUT=\
	"spline_Rk_integrals"


"$(INTDIR)\spline_Rk_integrals.obj"	"$(INTDIR)\spline_Rk_integrals.mod" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\spline_slater.f90
F90_MODOUT=\
	"spline_slater"


"$(INTDIR)\spline_slater.obj"	"$(INTDIR)\spline_slater.mod" : $(SOURCE) "$(INTDIR)"
	$(F90) $(F90_PROJ) $(SOURCE)


SOURCE=.\sum_amb.f90

"$(INTDIR)\sum_amb.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\tk.f90

"$(INTDIR)\tk.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\tk_moments.f90

"$(INTDIR)\tk_moments.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\tkc.f90

"$(INTDIR)\tkc.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\tky.f90

"$(INTDIR)\tky.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\vbsplvd.f90

"$(INTDIR)\vbsplvd.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\vinty.f90

"$(INTDIR)\vinty.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\vk.f90

"$(INTDIR)\vk.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\vk_moments.f90

"$(INTDIR)\vk_moments.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\vkc.f90

"$(INTDIR)\vkc.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\vky.f90

"$(INTDIR)\vky.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\wk.f90

"$(INTDIR)\wk.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\wk_moments.f90

"$(INTDIR)\wk_moments.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\wky.f90

"$(INTDIR)\wky.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\ykf.f90

"$(INTDIR)\ykf.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\yval.f90

"$(INTDIR)\yval.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\zeta_y.f90

"$(INTDIR)\zeta_y.obj" : $(SOURCE) "$(INTDIR)"



!ENDIF 

