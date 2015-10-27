# Microsoft Developer Studio Project File - Name="BS" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=BS - WIN32 RELEASE
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "BS.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "BS.mak" CFG="BS - WIN32 RELEASE"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "BS - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe
# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD F90 /compile_only /include:"Release/" /nologo /warn:argument_checking /warn:declarations /warn:nofileopt
# SUBTRACT F90 /check:bounds /check:overflow /traceback
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"C:\oleg\LIB\BS.lib"
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=copy release\*.mod C:\oleg\INCLUDE	copy release\*.mod .
# End Special Build Tool
# Begin Target

# Name "BS - Win32 Release"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\azl.f90
# End Source File
# Begin Source File

SOURCE=.\bav.f90
# End Source File
# Begin Source File

SOURCE=.\bcore.f90
DEP_F90_BCORE=\
	".\spline_atomic.mod"\
	".\spline_orbitals.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\bhl.f90
DEP_F90_BHL_F=\
	".\spline_atomic.mod"\
	".\spline_hl.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\bhwf.f90
DEP_F90_BHWF_=\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\bvalue2.f90
# End Source File
# Begin Source File

SOURCE=.\bvmv.f90
# End Source File
# Begin Source File

SOURCE=.\bxv.f90
# End Source File
# Begin Source File

SOURCE=.\bzk.f90
DEP_F90_BZK_F=\
	".\spline_grid.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	".\spline_slater.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\convol.f90
DEP_F90_CONVO=\
	".\spline_integrals.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\coulom.f90
DEP_F90_COULO=\
	".\spline_galerkin.mod"\
	".\spline_hl.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\define_grid.f90
DEP_F90_DEFIN=\
	".\spline_grid.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\define_spline.f90
DEP_F90_DEFINE=\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\density.f90
# End Source File
# Begin Source File

SOURCE=.\dinty.f90
DEP_F90_DINTY=\
	".\spline_grid.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\dvx.f90
# End Source File
# Begin Source File

SOURCE=.\facdyk.f90
DEP_F90_FACDY=\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\facdzk.f90
DEP_F90_FACDZ=\
	".\spline_galerkin.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\facsb.f90
DEP_F90_FACSB=\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\facsbl.f90
DEP_F90_FACSBL=\
	".\spline_galerkin.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\gauss.f90
# End Source File
# Begin Source File

SOURCE=.\grad.f90
DEP_F90_GRAD_=\
	".\spline_galerkin.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\hlc.f90
DEP_F90_HLC_F=\
	".\spline_atomic.mod"\
	".\spline_orbitals.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\hlm.f90
DEP_F90_HLM_F=\
	".\spline_atomic.mod"\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_hl.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\int_de.f90
DEP_F90_INT_D=\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\int_v.f90
DEP_F90_INT_V=\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\minty.f90
DEP_F90_MINTY=\
	".\spline_grid.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mk.f90
DEP_F90_MK_F9=\
	".\spline_integrals.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mk_moments.f90
DEP_F90_MK_MO=\
	".\spline_grid.mod"\
	".\spline_integrals.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mkc.f90
DEP_F90_MKC_F=\
	".\spline_atomic.mod"\
	".\spline_moments.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mky.f90
# End Source File
# Begin Source File

SOURCE=.\mmk_cell.f90
DEP_F90_MMK_C=\
	".\spline_atomic.mod"\
	".\spline_integrals.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mmk_diff.f90
DEP_F90_MMK_D=\
	".\spline_atomic.mod"\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_integrals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mnk_cell.f90
DEP_F90_MNK_C=\
	".\spline_atomic.mod"\
	".\spline_integrals.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mnk_diff.f90
DEP_F90_MNK_D=\
	".\spline_atomic.mod"\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_integrals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\moments.f90
DEP_F90_MOMEN=\
	".\spline_grid.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mqk_cell.f90
DEP_F90_MQK_C=\
	".\spline_atomic.mod"\
	".\spline_integrals.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mrk_cell.f90
DEP_F90_MRK_C=\
	".\spline_integrals.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mrk_diff.f90
DEP_F90_MRK_D=\
	".\spline_atomic.mod"\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_integrals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mrm.f90
DEP_F90_MRM_F=\
	".\spline_grid.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mtk_cell.f90
DEP_F90_MTK_C=\
	".\spline_atomic.mod"\
	".\spline_integrals.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mtk_diff.f90
DEP_F90_MTK_D=\
	".\spline_atomic.mod"\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_integrals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mvc.f90
DEP_F90_MVC_F=\
	".\spline_atomic.mod"\
	".\spline_grid.mod"\
	".\spline_hl.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mvcv.f90
DEP_F90_MVCV_=\
	".\spline_atomic.mod"\
	".\spline_grid.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mvk_cell.f90
DEP_F90_MVK_C=\
	".\spline_atomic.mod"\
	".\spline_integrals.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mvk_diff.f90
DEP_F90_MVK_D=\
	".\spline_atomic.mod"\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_integrals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mwk_cell.f90
DEP_F90_MWK_C=\
	".\spline_atomic.mod"\
	".\spline_integrals.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\nk.f90
DEP_F90_NK_F9=\
	".\spline_integrals.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\nk_moments.f90
DEP_F90_NK_MO=\
	".\spline_grid.mod"\
	".\spline_integrals.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\nkc.f90
DEP_F90_NKC_F=\
	".\spline_atomic.mod"\
	".\spline_moments.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\nky.f90
DEP_F90_NKY_F=\
	".\spline_atomic.mod"\
	".\spline_grid.mod"\
	".\spline_orbitals.mod"\
	".\spline_slater.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\qk.f90
DEP_F90_QK_F9=\
	".\spline_integrals.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\qk_moments.f90
DEP_F90_QK_MO=\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\qkc.f90
DEP_F90_QKC_F=\
	".\spline_atomic.mod"\
	".\spline_moments.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\qky.f90
# End Source File
# Begin Source File

SOURCE=.\quadr.f90
DEP_F90_QUADR=\
	".\spline_galerkin.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\r_bwfn.f90
DEP_F90_R_BWF=\
	".\spline_atomic.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\rk.f90
DEP_F90_RK_F9=\
	".\spline_integrals.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\rk_moments.f90
DEP_F90_RK_MO=\
	".\spline_atomic.mod"\
	".\spline_grid.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\rkc.f90
DEP_F90_RKC_F=\
	".\spline_moments.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\rky.f90
DEP_F90_RKY_F=\
	".\spline_atomic.mod"\
	".\spline_grid.mod"\
	".\spline_orbitals.mod"\
	".\spline_slater.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\splin3.f90
# End Source File
# Begin Source File

SOURCE=.\spline_atomic.f90
F90_MODOUT=\
	"spline_atomic"

# End Source File
# Begin Source File

SOURCE=.\spline_densities.f90
DEP_F90_SPLIN=\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
F90_MODOUT=\
	"spline_densities"

# End Source File
# Begin Source File

SOURCE=.\spline_galerkin.f90
DEP_F90_SPLINE=\
	".\spline_param.mod"\
	
F90_MODOUT=\
	"spline_galerkin"

# End Source File
# Begin Source File

SOURCE=.\spline_grid.f90
DEP_F90_SPLINE_=\
	".\spline_param.mod"\
	
# ADD F90 /assume:nosource_include
F90_MODOUT=\
	"spline_grid"

# End Source File
# Begin Source File

SOURCE=.\spline_hl.f90
DEP_F90_SPLINE_H=\
	".\spline_param.mod"\
	
F90_MODOUT=\
	"spline_hl"

# End Source File
# Begin Source File

SOURCE=.\spline_integrals.f90
DEP_F90_SPLINE_I=\
	".\spline_param.mod"\
	
F90_MODOUT=\
	"spline_integrals"

# End Source File
# Begin Source File

SOURCE=.\spline_moments.f90
DEP_F90_SPLINE_M=\
	".\spline_param.mod"\
	
F90_MODOUT=\
	"spline_moments"

# End Source File
# Begin Source File

SOURCE=.\spline_orbitals.f90
DEP_F90_SPLINE_O=\
	".\spline_param.mod"\
	
F90_MODOUT=\
	"spline_orbitals"

# End Source File
# Begin Source File

SOURCE=.\spline_param.f90
F90_MODOUT=\
	"spline_param"

# End Source File
# Begin Source File

SOURCE=.\spline_Rk_integrals.f90
DEP_F90_SPLINE_R=\
	".\spline_densities.mod"\
	".\spline_galerkin.mod"\
	".\spline_moments.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
F90_MODOUT=\
	"spline_Rk_integrals"

# End Source File
# Begin Source File

SOURCE=.\spline_slater.f90
DEP_F90_SPLINE_S=\
	".\spline_param.mod"\
	
F90_MODOUT=\
	"spline_slater"

# End Source File
# Begin Source File

SOURCE=.\sum_amb.f90
# End Source File
# Begin Source File

SOURCE=.\tk.f90
DEP_F90_TK_F9=\
	".\spline_integrals.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\tk_moments.f90
DEP_F90_TK_MO=\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\tkc.f90
DEP_F90_TKC_F=\
	".\spline_atomic.mod"\
	".\spline_moments.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\tky.f90
DEP_F90_TKY_F=\
	".\spline_atomic.mod"\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	".\spline_slater.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\vbsplvd.f90
DEP_F90_VBSPL=\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\vinty.f90
DEP_F90_VINTY=\
	".\spline_grid.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\vk.f90
DEP_F90_VK_F9=\
	".\spline_integrals.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\vk_moments.f90
DEP_F90_VK_MO=\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\vkc.f90
DEP_F90_VKC_F=\
	".\spline_atomic.mod"\
	".\spline_moments.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\vky.f90
DEP_F90_VKY_F=\
	".\spline_atomic.mod"\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	".\spline_slater.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\wk.f90
DEP_F90_WK_F9=\
	".\spline_integrals.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\wk_moments.f90
DEP_F90_WK_MO=\
	".\spline_galerkin.mod"\
	".\spline_grid.mod"\
	".\spline_moments.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\wky.f90
# End Source File
# Begin Source File

SOURCE=.\ykf.f90
DEP_F90_YKF_F=\
	".\spline_atomic.mod"\
	".\spline_grid.mod"\
	".\spline_orbitals.mod"\
	".\spline_param.mod"\
	".\spline_slater.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\yval.f90
DEP_F90_YVAL_=\
	".\spline_grid.mod"\
	".\spline_param.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\zeta_y.f90
DEP_F90_ZETA_=\
	".\spline_atomic.mod"\
	".\spline_orbitals.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# End Target
# End Project
