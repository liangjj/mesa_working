module blas95
  !================================================================================
  ! This module provides Fortran 95 wrappers to the legacy BLAS. It allows a
  ! Fortran programmer to use BLAS conveniently using array sections, whole
  ! arrays and optional arguments. The Level 2 and 3 interfaces here are BLAST
  ! compliant.
  !================================================================================

  implicit none
  private
  integer,parameter:: skind = kind(1.e0), dkind = kind(1.d0)

  ! ones
  real(skind),parameter:: s_one = real(1,kind=skind)
  real(dkind),parameter:: d_one = real(1,kind=dkind)
  complex(skind),parameter:: c_one = cmplx(1,kind=skind)
  complex(dkind),parameter:: z_one = cmplx(1,kind=dkind)
  ! zeros
  real(skind),parameter:: s_zero = real(0,kind=skind)
  real(dkind),parameter:: d_zero = real(0,kind=dkind)
  complex(skind),parameter:: c_zero = cmplx(0,kind=skind)
  complex(dkind),parameter:: z_zero = cmplx(0,kind=dkind)
  ! enable these if your compiler supports sizeof
  !integer,parameter:: ssize = sizeof(s_one)
  !integer,parameter:: dsize = sizeof(d_one)
  !integer,parameter:: csize = sizeof(c_one)
  !integer,parameter:: zsize = sizeof(z_one)
  ! sizeof is less common extension, so provide sizes manually
  integer,parameter:: ssize = 4 
  integer,parameter:: dsize = 8 
  integer,parameter:: csize = 2*ssize
  integer,parameter:: zsize = 2*dsize

  type,public:: blas_trans_type
     character:: c
  end type blas_trans_type
  type(blas_trans_type),parameter,public:: blas_trans = blas_trans_type("T")
  type(blas_trans_type),parameter,public:: blas_no_trans = blas_trans_type("N")
  type(blas_trans_type),parameter,public:: blas_conj_trans = blas_trans_type("H")

  type,public:: blas_uplo_type
     character:: c
  end type blas_uplo_type
  type(blas_uplo_type),parameter,public:: blas_lower = blas_uplo_type("L")
  type(blas_uplo_type),parameter,public:: blas_upper = blas_uplo_type("U")

  type,public:: blas_diag_type
     character:: c
  end type blas_diag_type
  type(blas_diag_type),parameter,public:: blas_unit_diag = blas_diag_type("U")
  type(blas_diag_type),parameter,public:: blas_non_unit_diag = blas_diag_type("N")

  type,public:: blas_side_type
     character:: c
  end type blas_side_type
  type(blas_side_type),parameter,public:: blas_left_side = blas_side_type("L")
  type(blas_side_type),parameter,public:: blas_right_side = blas_side_type("R")







  interface GEMM
     module procedure sGEMV_95
     module procedure dGEMV_95
     module procedure cGEMV_95
     module procedure zGEMV_95
  end interface


  public GEMM

  interface GBMV
     module procedure sGBMV_95
     module procedure dGBMV_95
     module procedure cGBMV_95
     module procedure zGBMV_95
  end interface


  public GBMV

  interface SYMM
     module procedure sSYMV_95
     module procedure dSYMV_95
  end interface


  public SYMM

  interface SBMV
     module procedure sSBMV_95
     module procedure dSBMV_95
  end interface


  public SBMV

  interface SPMV
     module procedure sSPMV_95
     module procedure dSPMV_95
  end interface


  public SPMV

  interface HEMM
     module procedure cHEMV_95
     module procedure zHEMV_95
  end interface


  public HEMM

  interface HBMV
     module procedure cHBMV_95
     module procedure zHBMV_95
  end interface


  public HBMV

  interface HPMV
     module procedure cHPMV_95
     module procedure zHPMV_95
  end interface


  public HPMV

  interface TRMM
     module procedure sTRMV_95
     module procedure dTRMV_95
     module procedure cTRMV_95
     module procedure zTRMV_95
  end interface


  public TRMM

  interface TBMV
     module procedure sTBMV_95
     module procedure dTBMV_95
     module procedure cTBMV_95
     module procedure zTBMV_95
  end interface


  public TBMV

  interface TPMV
     module procedure sTPMV_95
     module procedure dTPMV_95
     module procedure cTPMV_95
     module procedure zTPMV_95
  end interface


  public TPMV

  interface TRSM
     module procedure sTRSV_95
     module procedure dTRSV_95
     module procedure cTRSV_95
     module procedure zTRSV_95
  end interface


  public TRSM

  interface TBSV
     module procedure sTBSV_95
     module procedure dTBSV_95
     module procedure cTBSV_95
     module procedure zTBSV_95
  end interface


  public TBSV

  interface TPSV
     module procedure sTPSV_95
     module procedure dTPSV_95
     module procedure cTPSV_95
     module procedure zTPSV_95
  end interface


  public TPSV

  interface GEMM
     module procedure sGER_95
     module procedure dGER_95
  end interface

  ! GEMM already public
  interface GEMM
     module procedure cGERUC_95
     module procedure zGERUC_95
  end interface

  ! GEMM already public
  interface SYRK
     module procedure sSYR_95
     module procedure dSYR_95
  end interface


  public SYRK

  interface SPR
     module procedure sSPR_95
     module procedure dSPR_95
  end interface


  public SPR

  interface HERK
     module procedure cHER_95
     module procedure zHER_95
  end interface


  public HERK

  interface HPR
     module procedure cHPR_95
     module procedure zHPR_95
  end interface


  public HPR

  interface SYR2K
     module procedure sSYR2_95
     module procedure dSYR2_95
  end interface


  public SYR2K

  interface SPR2
     module procedure sSPR_95
     module procedure dSPR_95
  end interface


  public SPR2

  interface HER2K
     module procedure cHER2_95
     module procedure zHER2_95
  end interface


  public HER2K

  interface HPR2
     module procedure cHPR_95
     module procedure zHPR_95
  end interface


  public HPR2

  interface GEMM
     module procedure sGEMM_95
     module procedure dGEMM_95
     module procedure cGEMM_95
     module procedure zGEMM_95
  end interface

  ! GEMM already public
  interface SYMM
     module procedure sSYMM_95
     module procedure dSYMM_95
     module procedure cSYMM_95
     module procedure zSYMM_95
  end interface

  ! SYMM already public
  interface HEMM
     module procedure cHEMM_95
     module procedure zHEMM_95
  end interface

  ! HEMM already public
  interface TRMM
     module procedure sTRMM_95
     module procedure dTRMM_95
     module procedure cTRMM_95
     module procedure zTRMM_95
  end interface

  ! TRMM already public
  interface TRSM
     module procedure sTRSM_95
     module procedure dTRSM_95
     module procedure cTRSM_95
     module procedure zTRSM_95
  end interface

  ! TRSM already public
  interface SYRK
     module procedure sSYRK_95
     module procedure dSYRK_95
  end interface

  ! SYRK already public
  interface HERK
     module procedure cHERK_95
     module procedure zHERK_95
  end interface

  ! HERK already public
  interface SYR2K
     module procedure sSYR2K_95
     module procedure dSYR2K_95
  end interface

  ! SYR2K already public
  interface HER2K
     module procedure cHER2K_95
     module procedure zHER2K_95
  end interface

  ! HER2K already public

  public blas_OK_vector
  interface blas_OK_vector
     module procedure sblas_OK_vector
     module procedure dblas_OK_vector
     module procedure cblas_OK_vector
     module procedure zblas_OK_vector
  end interface
  public blas_OK_matrix
  interface blas_OK_matrix
     module procedure sblas_cont_vector
     module procedure dblas_cont_vector
     module procedure cblas_cont_vector
     module procedure zblas_cont_vector
     module procedure sblas_OK_matrix
     module procedure dblas_OK_matrix
     module procedure cblas_OK_matrix
     module procedure zblas_OK_matrix
  end interface

contains



  function sblas_OK_vector(vec) result(OK)
    real(skind),dimension(:):: vec
    logical:: OK
    OK = .true.
    if (size(vec) > 1) then
       OK = mod(loc(vec(2))-loc(vec(1)),ssize) == 0
    end if
  end function sblas_OK_vector

  function dblas_OK_vector(vec) result(OK)
    real(dkind),dimension(:):: vec
    logical:: OK
    OK = .true.
    if (size(vec) > 1) then
       OK = mod(loc(vec(2))-loc(vec(1)),dsize) == 0
    end if
  end function dblas_OK_vector

  function cblas_OK_vector(vec) result(OK)
    complex(skind),dimension(:):: vec
    logical:: OK
    OK = .true.
    if (size(vec) > 1) then
       OK = mod(loc(vec(2))-loc(vec(1)),csize) == 0
    end if
  end function cblas_OK_vector

  function zblas_OK_vector(vec) result(OK)
    complex(dkind),dimension(:):: vec
    logical:: OK
    OK = .true.
    if (size(vec) > 1) then
       OK = mod(loc(vec(2))-loc(vec(1)),zsize) == 0
    end if
  end function zblas_OK_vector



  function sblas_OK_matrix(mat) result(OK)
    real(skind),dimension(:,:):: mat
    logical:: OK
    OK = .true.
    if (size(mat,1) > 1 .and. size(mat,2) > 0) then
       OK = loc(mat(2,1))-loc(mat(1,1)) == ssize
    end if
    if (size(mat,2) > 1 .and. size(mat,1) > 0) then
       OK = OK .and. mod(loc(mat(1,2))-loc(mat(1,1)),ssize) == 0
    end if
  end function sblas_OK_matrix

  function dblas_OK_matrix(mat) result(OK)
    real(dkind),dimension(:,:):: mat
    logical:: OK
    OK = .true.
    if (size(mat,1) > 1 .and. size(mat,2) > 0) then
       OK = loc(mat(2,1))-loc(mat(1,1)) == dsize
    end if
    if (size(mat,2) > 1 .and. size(mat,1) > 0) then
       OK = OK .and. mod(loc(mat(1,2))-loc(mat(1,1)),dsize) == 0
    end if
  end function dblas_OK_matrix

  function cblas_OK_matrix(mat) result(OK)
    complex(skind),dimension(:,:):: mat
    logical:: OK
    OK = .true.
    if (size(mat,1) > 1 .and. size(mat,2) > 0) then
       OK = loc(mat(2,1))-loc(mat(1,1)) == csize
    end if
    if (size(mat,2) > 1 .and. size(mat,1) > 0) then
       OK = OK .and. mod(loc(mat(1,2))-loc(mat(1,1)),csize) == 0
    end if
  end function cblas_OK_matrix

  function zblas_OK_matrix(mat) result(OK)
    complex(dkind),dimension(:,:):: mat
    logical:: OK
    OK = .true.
    if (size(mat,1) > 1 .and. size(mat,2) > 0) then
       OK = loc(mat(2,1))-loc(mat(1,1)) == zsize
    end if
    if (size(mat,2) > 1 .and. size(mat,1) > 0) then
       OK = OK .and. mod(loc(mat(1,2))-loc(mat(1,1)),zsize) == 0
    end if
  end function zblas_OK_matrix



  function sblas_cont_vector(mat) result(OK)
    real(skind),dimension(:):: mat
    logical:: OK
    OK = .true.
    if (size(mat) > 1) OK = loc(mat(2))-loc(mat(1)) == ssize
  end function sblas_cont_vector

  function dblas_cont_vector(mat) result(OK)
    real(dkind),dimension(:):: mat
    logical:: OK
    OK = .true.
    if (size(mat) > 1) OK = loc(mat(2))-loc(mat(1)) == dsize
  end function dblas_cont_vector

  function cblas_cont_vector(mat) result(OK)
    complex(skind),dimension(:):: mat
    logical:: OK
    OK = .true.
    if (size(mat) > 1) OK = loc(mat(2))-loc(mat(1)) == csize
  end function cblas_cont_vector

  function zblas_cont_vector(mat) result(OK)
    complex(dkind),dimension(:):: mat
    logical:: OK
    OK = .true.
    if (size(mat) > 1) OK = loc(mat(2))-loc(mat(1)) == zsize
  end function zblas_cont_vector












  !================================================================================
  ! LEVEL 2 BLAS routines
  !================================================================================

  !========================================
  ! GEMV implementation



  subroutine sGEMV_95(a,b,c,transa,alpha,beta)
    real(skind),dimension(:,:),intent(in):: a
    real(skind),dimension(:),intent(in):: b
    real(skind),dimension(:),intent(inout):: c
    type(blas_trans_type),intent(in),optional:: transa
    character:: vtransa
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    real(skind),intent(in),optional:: beta
    real(skind):: vbeta
    integer:: lda,incb,incc
    integer:: m,n

    vtransa = "N"
    if (present(transa)) vtransa = transa%c

    valpha = s_one
    if (present(alpha)) valpha = alpha

    vbeta = s_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/ssize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/ssize
    incc = 1
    if (size(c) > 1) incc = (loc(c(2))-loc(c(1)))/ssize
    m = size(a,1)
    n = size(a,2)
    call sGEMV(vtransa,m,n,valpha,a(1,1),lda,b(1),incb,vbeta,c(1),incc)
  end subroutine sGEMV_95

  subroutine dGEMV_95(a,b,c,transa,alpha,beta)
    real(dkind),dimension(:,:),intent(in):: a
    real(dkind),dimension(:),intent(in):: b
    real(dkind),dimension(:),intent(inout):: c
    type(blas_trans_type),intent(in),optional:: transa
    character:: vtransa
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    real(dkind),intent(in),optional:: beta
    real(dkind):: vbeta
    integer:: lda,incb,incc
    integer:: m,n

    vtransa = "N"
    if (present(transa)) vtransa = transa%c

    valpha = d_one
    if (present(alpha)) valpha = alpha

    vbeta = d_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/dsize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/dsize
    incc = 1
    if (size(c) > 1) incc = (loc(c(2))-loc(c(1)))/dsize
    m = size(a,1)
    n = size(a,2)
    call dGEMV(vtransa,m,n,valpha,a(1,1),lda,b(1),incb,vbeta,c(1),incc)
  end subroutine dGEMV_95

  subroutine cGEMV_95(a,b,c,transa,alpha,beta)
    complex(skind),dimension(:,:),intent(in):: a
    complex(skind),dimension(:),intent(in):: b
    complex(skind),dimension(:),intent(inout):: c
    type(blas_trans_type),intent(in),optional:: transa
    character:: vtransa
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    complex(skind),intent(in),optional:: beta
    complex(skind):: vbeta
    integer:: lda,incb,incc
    integer:: m,n

    vtransa = "N"
    if (present(transa)) vtransa = transa%c

    valpha = c_one
    if (present(alpha)) valpha = alpha

    vbeta = c_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/csize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/csize
    incc = 1
    if (size(c) > 1) incc = (loc(c(2))-loc(c(1)))/csize
    m = size(a,1)
    n = size(a,2)
    call cGEMV(vtransa,m,n,valpha,a(1,1),lda,b(1),incb,vbeta,c(1),incc)
  end subroutine cGEMV_95

  subroutine zGEMV_95(a,b,c,transa,alpha,beta)
    complex(dkind),dimension(:,:),intent(in):: a
    complex(dkind),dimension(:),intent(in):: b
    complex(dkind),dimension(:),intent(inout):: c
    type(blas_trans_type),intent(in),optional:: transa
    character:: vtransa
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    complex(dkind),intent(in),optional:: beta
    complex(dkind):: vbeta
    integer:: lda,incb,incc
    integer:: m,n

    vtransa = "N"
    if (present(transa)) vtransa = transa%c

    valpha = z_one
    if (present(alpha)) valpha = alpha

    vbeta = z_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/zsize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/zsize
    incc = 1
    if (size(c) > 1) incc = (loc(c(2))-loc(c(1)))/zsize
    m = size(a,1)
    n = size(a,2)
    call zGEMV(vtransa,m,n,valpha,a(1,1),lda,b(1),incb,vbeta,c(1),incc)
  end subroutine zGEMV_95

  !========================================
  ! GBMV implementation


  subroutine sGBMV_95(a,m,kl,x,y,trans,alpha,beta)
    real(skind),dimension(:,:),intent(in):: a
    real(skind),dimension(:),intent(in):: x
    real(skind),dimension(:),intent(inout):: y
    integer,intent(in):: m,kl
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    real(skind),intent(in),optional:: beta
    real(skind):: vbeta
    integer:: n,ku,lda,incx,incy

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = s_one
    if (present(alpha)) valpha = alpha

    vbeta = s_zero
    if (present(beta)) vbeta = beta
    n = size(x)
    ku = size(a,1)-1-kl
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/ssize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/ssize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/ssize
    call sGBMV(vtrans,m,n,kl,ku,valpha,a(1,1),lda,x(1),incx,vbeta,y(1),incy)
  end subroutine sGBMV_95

  subroutine dGBMV_95(a,m,kl,x,y,trans,alpha,beta)
    real(dkind),dimension(:,:),intent(in):: a
    real(dkind),dimension(:),intent(in):: x
    real(dkind),dimension(:),intent(inout):: y
    integer,intent(in):: m,kl
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    real(dkind),intent(in),optional:: beta
    real(dkind):: vbeta
    integer:: n,ku,lda,incx,incy

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = d_one
    if (present(alpha)) valpha = alpha

    vbeta = d_zero
    if (present(beta)) vbeta = beta
    n = size(x)
    ku = size(a,1)-1-kl
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/dsize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/dsize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/dsize
    call dGBMV(vtrans,m,n,kl,ku,valpha,a(1,1),lda,x(1),incx,vbeta,y(1),incy)
  end subroutine dGBMV_95

  subroutine cGBMV_95(a,m,kl,x,y,trans,alpha,beta)
    complex(skind),dimension(:,:),intent(in):: a
    complex(skind),dimension(:),intent(in):: x
    complex(skind),dimension(:),intent(inout):: y
    integer,intent(in):: m,kl
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    complex(skind),intent(in),optional:: beta
    complex(skind):: vbeta
    integer:: n,ku,lda,incx,incy

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = c_one
    if (present(alpha)) valpha = alpha

    vbeta = c_zero
    if (present(beta)) vbeta = beta
    n = size(x)
    ku = size(a,1)-1-kl
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/csize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/csize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/csize
    call cGBMV(vtrans,m,n,kl,ku,valpha,a(1,1),lda,x(1),incx,vbeta,y(1),incy)
  end subroutine cGBMV_95

  subroutine zGBMV_95(a,m,kl,x,y,trans,alpha,beta)
    complex(dkind),dimension(:,:),intent(in):: a
    complex(dkind),dimension(:),intent(in):: x
    complex(dkind),dimension(:),intent(inout):: y
    integer,intent(in):: m,kl
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    complex(dkind),intent(in),optional:: beta
    complex(dkind):: vbeta
    integer:: n,ku,lda,incx,incy

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = z_one
    if (present(alpha)) valpha = alpha

    vbeta = z_zero
    if (present(beta)) vbeta = beta
    n = size(x)
    ku = size(a,1)-1-kl
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/zsize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/zsize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/zsize
    call zGBMV(vtrans,m,n,kl,ku,valpha,a(1,1),lda,x(1),incx,vbeta,y(1),incy)
  end subroutine zGBMV_95

  !========================================
  ! SYMV implementation



  subroutine sSYMV_95(a,b,c,uplo,alpha,beta)
    real(skind),dimension(:,:),intent(in):: a
    real(skind),dimension(:),intent(in):: b
    real(skind),dimension(:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    real(skind),intent(in),optional:: beta
    real(skind):: vbeta
    integer:: lda,incb,incc
    integer:: n

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = s_one
    if (present(alpha)) valpha = alpha

    vbeta = s_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/ssize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/ssize
    incc = 1
    if (size(c) > 1) incc = (loc(c(2))-loc(c(1)))/ssize
    n = size(a,2)
    call sSYMV(vuplo,n,valpha,a(1,1),lda,b(1),incb,vbeta,c(1),incc)
  end subroutine sSYMV_95

  subroutine dSYMV_95(a,b,c,uplo,alpha,beta)
    real(dkind),dimension(:,:),intent(in):: a
    real(dkind),dimension(:),intent(in):: b
    real(dkind),dimension(:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    real(dkind),intent(in),optional:: beta
    real(dkind):: vbeta
    integer:: lda,incb,incc
    integer:: n

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = d_one
    if (present(alpha)) valpha = alpha

    vbeta = d_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/dsize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/dsize
    incc = 1
    if (size(c) > 1) incc = (loc(c(2))-loc(c(1)))/dsize
    n = size(a,2)
    call dSYMV(vuplo,n,valpha,a(1,1),lda,b(1),incb,vbeta,c(1),incc)
  end subroutine dSYMV_95

  !========================================
  ! SBMV implementation


  subroutine sSBMV_95(a,x,y,uplo,alpha,beta)
    real(skind),dimension(:,:),intent(in):: a
    real(skind),dimension(:),intent(in):: x
    real(skind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    real(skind),intent(in),optional:: beta
    real(skind):: vbeta
    integer:: n,k,lda,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = s_one
    if (present(alpha)) valpha = alpha

    vbeta = s_zero
    if (present(beta)) vbeta = beta
    n = size(x)
    k = (size(a,1)-1)/2
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/ssize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/ssize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/ssize
    call sSBMV(vuplo,n,k,valpha,a(1,1),lda,x(1),incx,vbeta,y(1),incy)
  end subroutine sSBMV_95

  subroutine dSBMV_95(a,x,y,uplo,alpha,beta)
    real(dkind),dimension(:,:),intent(in):: a
    real(dkind),dimension(:),intent(in):: x
    real(dkind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    real(dkind),intent(in),optional:: beta
    real(dkind):: vbeta
    integer:: n,k,lda,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = d_one
    if (present(alpha)) valpha = alpha

    vbeta = d_zero
    if (present(beta)) vbeta = beta
    n = size(x)
    k = (size(a,1)-1)/2
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/dsize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/dsize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/dsize
    call dSBMV(vuplo,n,k,valpha,a(1,1),lda,x(1),incx,vbeta,y(1),incy)
  end subroutine dSBMV_95

  !========================================
  ! SPMV implementation


  subroutine sSPMV_95(ap,x,y,uplo,alpha,beta)
    real(skind),dimension(:),intent(in):: ap
    real(skind),dimension(:),intent(in):: x
    real(skind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    real(skind),intent(in),optional:: beta
    real(skind):: vbeta
    integer:: n,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = s_one
    if (present(alpha)) valpha = alpha

    vbeta = s_zero
    if (present(beta)) vbeta = beta
    n = size(x)
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/ssize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/ssize
    call sSPMV(vuplo,n,valpha,ap(1),x(1),incx,vbeta,y(1),incy)
  end subroutine sSPMV_95

  subroutine dSPMV_95(ap,x,y,uplo,alpha,beta)
    real(dkind),dimension(:),intent(in):: ap
    real(dkind),dimension(:),intent(in):: x
    real(dkind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    real(dkind),intent(in),optional:: beta
    real(dkind):: vbeta
    integer:: n,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = d_one
    if (present(alpha)) valpha = alpha

    vbeta = d_zero
    if (present(beta)) vbeta = beta
    n = size(x)
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/dsize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/dsize
    call dSPMV(vuplo,n,valpha,ap(1),x(1),incx,vbeta,y(1),incy)
  end subroutine dSPMV_95

  !========================================
  ! HEMV implementation



  subroutine cHEMV_95(a,b,c,uplo,alpha,beta)
    complex(skind),dimension(:,:),intent(in):: a
    complex(skind),dimension(:),intent(in):: b
    complex(skind),dimension(:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    complex(skind),intent(in),optional:: beta
    complex(skind):: vbeta
    integer:: lda,incb,incc
    integer:: n

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = c_one
    if (present(alpha)) valpha = alpha

    vbeta = c_zero
    if (present(beta)) vbeta = beta

    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/csize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/csize
    incc = 1
    if (size(c) > 1) incc = (loc(c(2))-loc(c(1)))/csize
    n = size(a,2)
    call cHEMV(vuplo,n,valpha,a(1,1),lda,b(1),incb,vbeta,c(1),incc)
  end subroutine cHEMV_95

  subroutine zHEMV_95(a,b,c,uplo,alpha,beta)
    complex(dkind),dimension(:,:),intent(in):: a
    complex(dkind),dimension(:),intent(in):: b
    complex(dkind),dimension(:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    complex(dkind),intent(in),optional:: beta
    complex(dkind):: vbeta
    integer:: lda,incb,incc
    integer:: n

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = z_one
    if (present(alpha)) valpha = alpha

    vbeta = z_zero
    if (present(beta)) vbeta = beta

    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/zsize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/zsize
    incc = 1
    if (size(c) > 1) incc = (loc(c(2))-loc(c(1)))/zsize
    n = size(a,2)
    call zHEMV(vuplo,n,valpha,a(1,1),lda,b(1),incb,vbeta,c(1),incc)
  end subroutine zHEMV_95

  !========================================
  ! HBMV implementation


  subroutine cHBMV_95(a,x,y,uplo,alpha,beta)
    complex(skind),dimension(:,:),intent(in):: a
    complex(skind),dimension(:),intent(in):: x
    complex(skind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    complex(skind),intent(in),optional:: beta
    complex(skind):: vbeta
    integer:: n,k,lda,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = c_one
    if (present(alpha)) valpha = alpha

    vbeta = c_zero
    if (present(beta)) vbeta = beta
    n = size(x)
    k = (size(a,1)-1)/2
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/csize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/csize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/csize
    call cHBMV(vuplo,n,k,valpha,a(1,1),lda,x(1),incx,vbeta,y(1),incy)
  end subroutine cHBMV_95

  subroutine zHBMV_95(a,x,y,uplo,alpha,beta)
    complex(dkind),dimension(:,:),intent(in):: a
    complex(dkind),dimension(:),intent(in):: x
    complex(dkind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    complex(dkind),intent(in),optional:: beta
    complex(dkind):: vbeta
    integer:: n,k,lda,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = z_one
    if (present(alpha)) valpha = alpha

    vbeta = z_zero
    if (present(beta)) vbeta = beta
    n = size(x)
    k = (size(a,1)-1)/2
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/zsize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/zsize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/zsize
    call zHBMV(vuplo,n,k,valpha,a(1,1),lda,x(1),incx,vbeta,y(1),incy)
  end subroutine zHBMV_95

  !========================================
  ! HPMV implementation


  subroutine cHPMV_95(ap,x,y,uplo,alpha,beta)
    complex(skind),dimension(:),intent(in):: ap
    complex(skind),dimension(:),intent(in):: x
    complex(skind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    complex(skind),intent(in),optional:: beta
    complex(skind):: vbeta
    integer:: n,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = c_one
    if (present(alpha)) valpha = alpha

    vbeta = c_zero
    if (present(beta)) vbeta = beta
    n = size(x)
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/csize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/csize
    call cHPMV(vuplo,n,valpha,ap(1),x(1),incx,vbeta,y(1),incy)
  end subroutine cHPMV_95

  subroutine zHPMV_95(ap,x,y,uplo,alpha,beta)
    complex(dkind),dimension(:),intent(in):: ap
    complex(dkind),dimension(:),intent(in):: x
    complex(dkind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    complex(dkind),intent(in),optional:: beta
    complex(dkind):: vbeta
    integer:: n,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = z_one
    if (present(alpha)) valpha = alpha

    vbeta = z_zero
    if (present(beta)) vbeta = beta
    n = size(x)
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/zsize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/zsize
    call zHPMV(vuplo,n,valpha,ap(1),x(1),incx,vbeta,y(1),incy)
  end subroutine zHPMV_95

  !========================================
  ! TRMV implementation



  subroutine sTRMV_95(t,b,side,uplo,transt,diag)
    real(skind),dimension(:,:),intent(in):: t
    real(skind),dimension(:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: ldt,incb
    integer:: n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/ssize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/ssize
    n = size(t,1)
    call sTRMV(vuplo,vtranst,vdiag,n,t(1,1),ldt,b(1),incb)
  end subroutine sTRMV_95

  subroutine dTRMV_95(t,b,side,uplo,transt,diag)
    real(dkind),dimension(:,:),intent(in):: t
    real(dkind),dimension(:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: ldt,incb
    integer:: n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/dsize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/dsize
    n = size(t,1)
    call dTRMV(vuplo,vtranst,vdiag,n,t(1,1),ldt,b(1),incb)
  end subroutine dTRMV_95

  subroutine cTRMV_95(t,b,side,uplo,transt,diag)
    complex(skind),dimension(:,:),intent(in):: t
    complex(skind),dimension(:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: ldt,incb
    integer:: n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/csize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/csize
    n = size(t,1)
    call cTRMV(vuplo,vtranst,vdiag,n,t(1,1),ldt,b(1),incb)
  end subroutine cTRMV_95

  subroutine zTRMV_95(t,b,side,uplo,transt,diag)
    complex(dkind),dimension(:,:),intent(in):: t
    complex(dkind),dimension(:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: ldt,incb
    integer:: n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/zsize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/zsize
    n = size(t,1)
    call zTRMV(vuplo,vtranst,vdiag,n,t(1,1),ldt,b(1),incb)
  end subroutine zTRMV_95

  !========================================
  ! TBMV implementation


  subroutine sTBMV_95(a,x,y,uplo,trans,diag)
    real(skind),dimension(:,:),intent(in):: a
    real(skind),dimension(:),intent(in):: x
    real(skind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,lda,incx

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    k = size(a,1)-1
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/ssize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/ssize
    call sTBMV(vuplo,vtrans,vdiag,n,k,a(1,1),lda,x(1),incx)
  end subroutine sTBMV_95

  subroutine dTBMV_95(a,x,y,uplo,trans,diag)
    real(dkind),dimension(:,:),intent(in):: a
    real(dkind),dimension(:),intent(in):: x
    real(dkind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,lda,incx

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    k = size(a,1)-1
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/dsize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/dsize
    call dTBMV(vuplo,vtrans,vdiag,n,k,a(1,1),lda,x(1),incx)
  end subroutine dTBMV_95

  subroutine cTBMV_95(a,x,y,uplo,trans,diag)
    complex(skind),dimension(:,:),intent(in):: a
    complex(skind),dimension(:),intent(in):: x
    complex(skind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,lda,incx

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    k = size(a,1)-1
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/csize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/csize
    call cTBMV(vuplo,vtrans,vdiag,n,k,a(1,1),lda,x(1),incx)
  end subroutine cTBMV_95

  subroutine zTBMV_95(a,x,y,uplo,trans,diag)
    complex(dkind),dimension(:,:),intent(in):: a
    complex(dkind),dimension(:),intent(in):: x
    complex(dkind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,lda,incx

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    k = size(a,1)-1
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/zsize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/zsize
    call zTBMV(vuplo,vtrans,vdiag,n,k,a(1,1),lda,x(1),incx)
  end subroutine zTBMV_95

  !========================================
  ! TPMV implementation


  subroutine sTPMV_95(ap,x,y,uplo,trans,diag)
    real(skind),dimension(:),intent(in):: ap
    real(skind),dimension(:),intent(in):: x
    real(skind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/ssize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/ssize
    call sTPMV(vuplo,vtrans,vdiag,n,ap(1),x(1),incx,y(1),incy)
  end subroutine sTPMV_95

  subroutine dTPMV_95(ap,x,y,uplo,trans,diag)
    real(dkind),dimension(:),intent(in):: ap
    real(dkind),dimension(:),intent(in):: x
    real(dkind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/dsize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/dsize
    call dTPMV(vuplo,vtrans,vdiag,n,ap(1),x(1),incx,y(1),incy)
  end subroutine dTPMV_95

  subroutine cTPMV_95(ap,x,y,uplo,trans,diag)
    complex(skind),dimension(:),intent(in):: ap
    complex(skind),dimension(:),intent(in):: x
    complex(skind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/csize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/csize
    call cTPMV(vuplo,vtrans,vdiag,n,ap(1),x(1),incx,y(1),incy)
  end subroutine cTPMV_95

  subroutine zTPMV_95(ap,x,y,uplo,trans,diag)
    complex(dkind),dimension(:),intent(in):: ap
    complex(dkind),dimension(:),intent(in):: x
    complex(dkind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/zsize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/zsize
    call zTPMV(vuplo,vtrans,vdiag,n,ap(1),x(1),incx,y(1),incy)
  end subroutine zTPMV_95

  !========================================
  ! TRSV implementation



  subroutine sTRSV_95(t,b,side,uplo,transt,diag)
    real(skind),dimension(:,:),intent(in):: t
    real(skind),dimension(:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: ldt,incb
    integer:: n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/ssize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/ssize
    n = size(t,1)
    call sTRSV(vuplo,vtranst,vdiag,n,t(1,1),ldt,b(1),incb)
  end subroutine sTRSV_95

  subroutine dTRSV_95(t,b,side,uplo,transt,diag)
    real(dkind),dimension(:,:),intent(in):: t
    real(dkind),dimension(:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: ldt,incb
    integer:: n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/dsize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/dsize
    n = size(t,1)
    call dTRSV(vuplo,vtranst,vdiag,n,t(1,1),ldt,b(1),incb)
  end subroutine dTRSV_95

  subroutine cTRSV_95(t,b,side,uplo,transt,diag)
    complex(skind),dimension(:,:),intent(in):: t
    complex(skind),dimension(:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: ldt,incb
    integer:: n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/csize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/csize
    n = size(t,1)
    call cTRSV(vuplo,vtranst,vdiag,n,t(1,1),ldt,b(1),incb)
  end subroutine cTRSV_95

  subroutine zTRSV_95(t,b,side,uplo,transt,diag)
    complex(dkind),dimension(:,:),intent(in):: t
    complex(dkind),dimension(:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: ldt,incb
    integer:: n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/zsize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/zsize
    n = size(t,1)
    call zTRSV(vuplo,vtranst,vdiag,n,t(1,1),ldt,b(1),incb)
  end subroutine zTRSV_95

  !========================================
  ! TBSV implementation


  subroutine sTBSV_95(a,x,y,uplo,trans,diag)
    real(skind),dimension(:,:),intent(in):: a
    real(skind),dimension(:),intent(in):: x
    real(skind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,lda,incx

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    k = size(a,1)-1
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/ssize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/ssize
    call sTBSV(vuplo,vtrans,vdiag,n,k,a(1,1),lda,x(1),incx)
  end subroutine sTBSV_95

  subroutine dTBSV_95(a,x,y,uplo,trans,diag)
    real(dkind),dimension(:,:),intent(in):: a
    real(dkind),dimension(:),intent(in):: x
    real(dkind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,lda,incx

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    k = size(a,1)-1
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/dsize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/dsize
    call dTBSV(vuplo,vtrans,vdiag,n,k,a(1,1),lda,x(1),incx)
  end subroutine dTBSV_95

  subroutine cTBSV_95(a,x,y,uplo,trans,diag)
    complex(skind),dimension(:,:),intent(in):: a
    complex(skind),dimension(:),intent(in):: x
    complex(skind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,lda,incx

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    k = size(a,1)-1
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/csize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/csize
    call cTBSV(vuplo,vtrans,vdiag,n,k,a(1,1),lda,x(1),incx)
  end subroutine cTBSV_95

  subroutine zTBSV_95(a,x,y,uplo,trans,diag)
    complex(dkind),dimension(:,:),intent(in):: a
    complex(dkind),dimension(:),intent(in):: x
    complex(dkind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,lda,incx

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    k = size(a,1)-1
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/zsize
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/zsize
    call zTBSV(vuplo,vtrans,vdiag,n,k,a(1,1),lda,x(1),incx)
  end subroutine zTBSV_95

  !========================================
  ! TPSV implementation


  subroutine sTPSV_95(ap,x,y,uplo,trans,diag)
    real(skind),dimension(:),intent(in):: ap
    real(skind),dimension(:),intent(in):: x
    real(skind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/ssize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/ssize
    call sTPSV(vuplo,vtrans,vdiag,n,ap(1),x(1),incx,y(1),incy)
  end subroutine sTPSV_95

  subroutine dTPSV_95(ap,x,y,uplo,trans,diag)
    real(dkind),dimension(:),intent(in):: ap
    real(dkind),dimension(:),intent(in):: x
    real(dkind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/dsize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/dsize
    call dTPSV(vuplo,vtrans,vdiag,n,ap(1),x(1),incx,y(1),incy)
  end subroutine dTPSV_95

  subroutine cTPSV_95(ap,x,y,uplo,trans,diag)
    complex(skind),dimension(:),intent(in):: ap
    complex(skind),dimension(:),intent(in):: x
    complex(skind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/csize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/csize
    call cTPSV(vuplo,vtrans,vdiag,n,ap(1),x(1),incx,y(1),incy)
  end subroutine cTPSV_95

  subroutine zTPSV_95(ap,x,y,uplo,trans,diag)
    complex(dkind),dimension(:),intent(in):: ap
    complex(dkind),dimension(:),intent(in):: x
    complex(dkind),dimension(:),intent(inout):: y
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    integer:: n,k,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c
    n = size(x)
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/zsize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/zsize
    call zTPSV(vuplo,vtrans,vdiag,n,ap(1),x(1),incx,y(1),incy)
  end subroutine zTPSV_95

  !========================================
  ! GER implementation



  subroutine sGER_95(a,b,c,alpha)
    real(skind),dimension(:),intent(in):: a
    real(skind),dimension(:),intent(in):: b
    real(skind),dimension(:,:),intent(inout):: c
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    integer:: ldc,inca,incb
    integer:: m,n

    valpha = s_one
    if (present(alpha)) valpha = alpha
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/ssize
    inca = 1
    if (size(a) > 1) inca = (loc(a(2))-loc(a(1)))/ssize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/ssize
    m = size(c,1)
    n = size(c,2)
    call sGER(m,n,valpha,a(1),inca,b(1),incb,c(1,1),ldc)
  end subroutine sGER_95

  subroutine dGER_95(a,b,c,alpha)
    real(dkind),dimension(:),intent(in):: a
    real(dkind),dimension(:),intent(in):: b
    real(dkind),dimension(:,:),intent(inout):: c
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    integer:: ldc,inca,incb
    integer:: m,n

    valpha = d_one
    if (present(alpha)) valpha = alpha
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/dsize
    inca = 1
    if (size(a) > 1) inca = (loc(a(2))-loc(a(1)))/dsize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/dsize
    m = size(c,1)
    n = size(c,2)
    call dGER(m,n,valpha,a(1),inca,b(1),incb,c(1,1),ldc)
  end subroutine dGER_95




  subroutine cGERUC_95(a,b,c,transb,alpha)
    complex(skind),dimension(:),intent(in):: a
    complex(skind),dimension(:),intent(in):: b
    complex(skind),dimension(:,:),intent(inout):: c
    type(blas_trans_type),intent(in),optional:: transb
    character:: vtransb
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    integer:: ldc,inca,incb
    integer:: m,n

    vtransb = "T"
    if (present(transb)) vtransb = transb%c

    valpha = c_one
    if (present(alpha)) valpha = alpha
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/csize
    inca = 1
    if (size(a) > 1) inca = (loc(a(2))-loc(a(1)))/csize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/csize
    m = size(c,1)
    n = size(c,2)
    if (vtransb == "H") then
       call cGERC(m,n,valpha,a(1),inca,b(1),incb,c(1,1),ldc)
    else 
       call cGERU(m,n,valpha,a(1),inca,b(1),incb,c(1,1),ldc)
    end if
  end subroutine cGERUC_95

  subroutine zGERUC_95(a,b,c,transb,alpha)
    complex(dkind),dimension(:),intent(in):: a
    complex(dkind),dimension(:),intent(in):: b
    complex(dkind),dimension(:,:),intent(inout):: c
    type(blas_trans_type),intent(in),optional:: transb
    character:: vtransb
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    integer:: ldc,inca,incb
    integer:: m,n

    vtransb = "T"
    if (present(transb)) vtransb = transb%c

    valpha = z_one
    if (present(alpha)) valpha = alpha
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/zsize
    inca = 1
    if (size(a) > 1) inca = (loc(a(2))-loc(a(1)))/zsize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/zsize
    m = size(c,1)
    n = size(c,2)
    if (vtransb == "H") then
       call zGERC(m,n,valpha,a(1),inca,b(1),incb,c(1,1),ldc)
    else 
       call zGERU(m,n,valpha,a(1),inca,b(1),incb,c(1,1),ldc)
    end if
  end subroutine zGERUC_95

  !========================================
  ! SYR implementation


  subroutine sSYR_95(a,c,uplo,alpha)
    real(skind),dimension(:),intent(in):: a
    real(skind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    integer:: n,ldc,inca

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = s_one
    if (present(alpha)) valpha = alpha
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/ssize
    inca = 1
    if (size(a) > 1) inca = (loc(a(2))-loc(a(1)))/ssize
    n = size(a)
    call sSYR(vuplo,n,valpha,a(1),inca,c(1,1),ldc)
  end subroutine sSYR_95

  subroutine dSYR_95(a,c,uplo,alpha)
    real(dkind),dimension(:),intent(in):: a
    real(dkind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    integer:: n,ldc,inca

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = d_one
    if (present(alpha)) valpha = alpha
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/dsize
    inca = 1
    if (size(a) > 1) inca = (loc(a(2))-loc(a(1)))/dsize
    n = size(a)
    call dSYR(vuplo,n,valpha,a(1),inca,c(1,1),ldc)
  end subroutine dSYR_95

  !========================================
  ! SPR implementation


  subroutine sSPR_95(x,ap,uplo,alpha)
    real(skind),dimension(:),intent(in):: x
    real(skind),dimension(:),intent(inout):: ap
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    integer:: n,incx

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = s_one
    if (present(alpha)) valpha = alpha
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/ssize
    n = size(x)
    call sSPR(vuplo,n,valpha,x(1),incx,ap(1))
  end subroutine sSPR_95

  subroutine dSPR_95(x,ap,uplo,alpha)
    real(dkind),dimension(:),intent(in):: x
    real(dkind),dimension(:),intent(inout):: ap
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    integer:: n,incx

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = d_one
    if (present(alpha)) valpha = alpha
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/dsize
    n = size(x)
    call dSPR(vuplo,n,valpha,x(1),incx,ap(1))
  end subroutine dSPR_95

  !========================================
  ! HER implementation


  subroutine cHER_95(a,c,uplo,alpha)
    complex(skind),dimension(:),intent(in):: a
    complex(skind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    integer:: n,ldc,inca

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = c_one
    if (present(alpha)) valpha = alpha
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/csize
    inca = 1
    if (size(a) > 1) inca = (loc(a(2))-loc(a(1)))/csize
    n = size(a)
    call cHER(vuplo,n,valpha,a(1),inca,c(1,1),ldc)
  end subroutine cHER_95

  subroutine zHER_95(a,c,uplo,alpha)
    complex(dkind),dimension(:),intent(in):: a
    complex(dkind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    integer:: n,ldc,inca

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = z_one
    if (present(alpha)) valpha = alpha
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/zsize
    inca = 1
    if (size(a) > 1) inca = (loc(a(2))-loc(a(1)))/zsize
    n = size(a)
    call zHER(vuplo,n,valpha,a(1),inca,c(1,1),ldc)
  end subroutine zHER_95

  !========================================
  ! HPR implementation


  subroutine cHPR_95(x,ap,uplo,alpha)
    complex(skind),dimension(:),intent(in):: x
    complex(skind),dimension(:),intent(inout):: ap
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    integer:: n,incx

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = c_one
    if (present(alpha)) valpha = alpha
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/csize
    n = size(x)
    call cHPR(vuplo,n,valpha,x(1),incx,ap(1))
  end subroutine cHPR_95

  subroutine zHPR_95(x,ap,uplo,alpha)
    complex(dkind),dimension(:),intent(in):: x
    complex(dkind),dimension(:),intent(inout):: ap
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    integer:: n,incx

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = z_one
    if (present(alpha)) valpha = alpha
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/zsize
    n = size(x)
    call zHPR(vuplo,n,valpha,x(1),incx,ap(1))
  end subroutine zHPR_95

  !========================================
  ! SYR2 implementation


  subroutine sSYR2_95(a,b,c,uplo,alpha)
    real(skind),dimension(:),intent(in):: a,b
    real(skind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    integer:: n,inca,incb,ldc

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = s_one
    if (present(alpha)) valpha = alpha
    inca = 1
    if (size(a) > 1) inca = (loc(a(2))-loc(a(1)))/ssize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/ssize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/ssize
    n = size(c,1)
    call sSYR2(vuplo,n,valpha,a(1),inca,b(1),incb,c(1,1),ldc)
  end subroutine sSYR2_95

  subroutine dSYR2_95(a,b,c,uplo,alpha)
    real(dkind),dimension(:),intent(in):: a,b
    real(dkind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    integer:: n,inca,incb,ldc

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = d_one
    if (present(alpha)) valpha = alpha
    inca = 1
    if (size(a) > 1) inca = (loc(a(2))-loc(a(1)))/dsize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/dsize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/dsize
    n = size(c,1)
    call dSYR2(vuplo,n,valpha,a(1),inca,b(1),incb,c(1,1),ldc)
  end subroutine dSYR2_95

  !========================================
  ! SPR2 implementation


  subroutine sSPR2_95(x,y,ap,uplo,alpha)
    real(skind),dimension(:),intent(in):: x,y
    real(skind),dimension(:),intent(inout):: ap
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    integer:: n,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = s_one
    if (present(alpha)) valpha = alpha
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/ssize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/ssize
    n = size(x)
    call sSPR2(vuplo,n,valpha,x(1),incx,y(1),incy,ap(1))
  end subroutine sSPR2_95

  subroutine dSPR2_95(x,y,ap,uplo,alpha)
    real(dkind),dimension(:),intent(in):: x,y
    real(dkind),dimension(:),intent(inout):: ap
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    integer:: n,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = d_one
    if (present(alpha)) valpha = alpha
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/dsize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/dsize
    n = size(x)
    call dSPR2(vuplo,n,valpha,x(1),incx,y(1),incy,ap(1))
  end subroutine dSPR2_95

  !========================================
  ! HER2 implementation


  subroutine cHER2_95(a,b,c,uplo,alpha)
    complex(skind),dimension(:),intent(in):: a,b
    complex(skind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    integer:: n,inca,incb,ldc

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = c_one
    if (present(alpha)) valpha = alpha
    inca = 1
    if (size(a) > 1) inca = (loc(a(2))-loc(a(1)))/csize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/csize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/csize
    n = size(c,1)
    call cHER2(vuplo,n,valpha,a(1),inca,b(1),incb,c(1,1),ldc)
  end subroutine cHER2_95

  subroutine zHER2_95(a,b,c,uplo,alpha)
    complex(dkind),dimension(:),intent(in):: a,b
    complex(dkind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    integer:: n,inca,incb,ldc

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = z_one
    if (present(alpha)) valpha = alpha
    inca = 1
    if (size(a) > 1) inca = (loc(a(2))-loc(a(1)))/zsize
    incb = 1
    if (size(b) > 1) incb = (loc(b(2))-loc(b(1)))/zsize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/zsize
    n = size(c,1)
    call zHER2(vuplo,n,valpha,a(1),inca,b(1),incb,c(1,1),ldc)
  end subroutine zHER2_95

  !========================================
  ! HPR2 implementation


  subroutine cHPR2_95(x,y,ap,uplo,alpha)
    complex(skind),dimension(:),intent(in):: x,y
    complex(skind),dimension(:),intent(inout):: ap
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    integer:: n,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = c_one
    if (present(alpha)) valpha = alpha
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/csize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/csize
    n = size(x)
    call cHPR2(vuplo,n,valpha,x(1),incx,y(1),incy,ap(1))
  end subroutine cHPR2_95

  subroutine zHPR2_95(x,y,ap,uplo,alpha)
    complex(dkind),dimension(:),intent(in):: x,y
    complex(dkind),dimension(:),intent(inout):: ap
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    integer:: n,incx,incy

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = z_one
    if (present(alpha)) valpha = alpha
    incx = 1
    if (size(x) > 1) incx = (loc(x(2))-loc(x(1)))/zsize
    incy = 1
    if (size(y) > 1) incy = (loc(y(2))-loc(y(1)))/zsize
    n = size(x)
    call zHPR2(vuplo,n,valpha,x(1),incx,y(1),incy,ap(1))
  end subroutine zHPR2_95


  !================================================================================
  ! LEVEL 3 BLAS routines
  !================================================================================

  !========================================
  ! GEMM implementation



  subroutine sGEMM_95(a,b,c,transa,transb,alpha,beta)
    real(skind),dimension(:,:),intent(in):: a
    real(skind),dimension(:,:),intent(in):: b
    real(skind),dimension(:,:),intent(inout):: c
    type(blas_trans_type),intent(in),optional:: transa
    character:: vtransa
    type(blas_trans_type),intent(in),optional:: transb
    character:: vtransb
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    real(skind),intent(in),optional:: beta
    real(skind):: vbeta
    integer:: lda,ldb,ldc
    integer:: m,n,k

    vtransa = "N"
    if (present(transa)) vtransa = transa%c

    vtransb = "N"
    if (present(transb)) vtransb = transb%c

    valpha = s_one
    if (present(alpha)) valpha = alpha

    vbeta = s_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/ssize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/ssize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/ssize
    m = size(c,1)
    n = size(c,2)
    k = merge(size(a,1),size(a,2),vtransa == "T")
    call sGEMM(vtransa,vtransb,m,n,k,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine sGEMM_95

  subroutine dGEMM_95(a,b,c,transa,transb,alpha,beta)
    real(dkind),dimension(:,:),intent(in):: a
    real(dkind),dimension(:,:),intent(in):: b
    real(dkind),dimension(:,:),intent(inout):: c
    type(blas_trans_type),intent(in),optional:: transa
    character:: vtransa
    type(blas_trans_type),intent(in),optional:: transb
    character:: vtransb
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    real(dkind),intent(in),optional:: beta
    real(dkind):: vbeta
    integer:: lda,ldb,ldc
    integer:: m,n,k

    vtransa = "N"
    if (present(transa)) vtransa = transa%c

    vtransb = "N"
    if (present(transb)) vtransb = transb%c

    valpha = d_one
    if (present(alpha)) valpha = alpha

    vbeta = d_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/dsize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/dsize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/dsize
    m = size(c,1)
    n = size(c,2)
    k = merge(size(a,1),size(a,2),vtransa == "T")
    call dGEMM(vtransa,vtransb,m,n,k,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine dGEMM_95

  subroutine cGEMM_95(a,b,c,transa,transb,alpha,beta)
    complex(skind),dimension(:,:),intent(in):: a
    complex(skind),dimension(:,:),intent(in):: b
    complex(skind),dimension(:,:),intent(inout):: c
    type(blas_trans_type),intent(in),optional:: transa
    character:: vtransa
    type(blas_trans_type),intent(in),optional:: transb
    character:: vtransb
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    complex(skind),intent(in),optional:: beta
    complex(skind):: vbeta
    integer:: lda,ldb,ldc
    integer:: m,n,k

    vtransa = "N"
    if (present(transa)) vtransa = transa%c

    vtransb = "N"
    if (present(transb)) vtransb = transb%c

    valpha = c_one
    if (present(alpha)) valpha = alpha

    vbeta = c_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/csize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/csize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/csize
    m = size(c,1)
    n = size(c,2)
    k = merge(size(a,1),size(a,2),vtransa == "T")
    call cGEMM(vtransa,vtransb,m,n,k,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine cGEMM_95

  subroutine zGEMM_95(a,b,c,transa,transb,alpha,beta)
    complex(dkind),dimension(:,:),intent(in):: a
    complex(dkind),dimension(:,:),intent(in):: b
    complex(dkind),dimension(:,:),intent(inout):: c
    type(blas_trans_type),intent(in),optional:: transa
    character:: vtransa
    type(blas_trans_type),intent(in),optional:: transb
    character:: vtransb
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    complex(dkind),intent(in),optional:: beta
    complex(dkind):: vbeta
    integer:: lda,ldb,ldc
    integer:: m,n,k

    vtransa = "N"
    if (present(transa)) vtransa = transa%c

    vtransb = "N"
    if (present(transb)) vtransb = transb%c

    valpha = z_one
    if (present(alpha)) valpha = alpha

    vbeta = z_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/zsize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/zsize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/zsize
    m = size(c,1)
    n = size(c,2)
    k = merge(size(a,1),size(a,2),vtransa == "T")
    call zGEMM(vtransa,vtransb,m,n,k,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine zGEMM_95

  !========================================
  ! SYMM implementation


  subroutine sSYMM_95(a,b,c,side,uplo,alpha,beta)
    real(skind),dimension(:,:),intent(in):: a
    real(skind),dimension(:,:),intent(in):: b
    real(skind),dimension(:,:),intent(inout):: c
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    real(skind),intent(in),optional:: beta
    real(skind):: vbeta
    integer:: lda,ldb,ldc
    integer:: m,n,k

    vside = "L"
    if (present(side)) vside = side%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = s_one
    if (present(alpha)) valpha = alpha

    vbeta = s_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/ssize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/ssize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/ssize
    m = size(c,1)
    n = size(c,2)
    call sSYMM(vside,vuplo,m,n,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine sSYMM_95

  subroutine dSYMM_95(a,b,c,side,uplo,alpha,beta)
    real(dkind),dimension(:,:),intent(in):: a
    real(dkind),dimension(:,:),intent(in):: b
    real(dkind),dimension(:,:),intent(inout):: c
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    real(dkind),intent(in),optional:: beta
    real(dkind):: vbeta
    integer:: lda,ldb,ldc
    integer:: m,n,k

    vside = "L"
    if (present(side)) vside = side%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = d_one
    if (present(alpha)) valpha = alpha

    vbeta = d_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/dsize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/dsize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/dsize
    m = size(c,1)
    n = size(c,2)
    call dSYMM(vside,vuplo,m,n,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine dSYMM_95

  subroutine cSYMM_95(a,b,c,side,uplo,alpha,beta)
    complex(skind),dimension(:,:),intent(in):: a
    complex(skind),dimension(:,:),intent(in):: b
    complex(skind),dimension(:,:),intent(inout):: c
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    complex(skind),intent(in),optional:: beta
    complex(skind):: vbeta
    integer:: lda,ldb,ldc
    integer:: m,n,k

    vside = "L"
    if (present(side)) vside = side%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = c_one
    if (present(alpha)) valpha = alpha

    vbeta = c_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/csize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/csize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/csize
    m = size(c,1)
    n = size(c,2)
    call cSYMM(vside,vuplo,m,n,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine cSYMM_95

  subroutine zSYMM_95(a,b,c,side,uplo,alpha,beta)
    complex(dkind),dimension(:,:),intent(in):: a
    complex(dkind),dimension(:,:),intent(in):: b
    complex(dkind),dimension(:,:),intent(inout):: c
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    complex(dkind),intent(in),optional:: beta
    complex(dkind):: vbeta
    integer:: lda,ldb,ldc
    integer:: m,n,k

    vside = "L"
    if (present(side)) vside = side%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = z_one
    if (present(alpha)) valpha = alpha

    vbeta = z_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/zsize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/zsize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/zsize
    m = size(c,1)
    n = size(c,2)
    call zSYMM(vside,vuplo,m,n,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine zSYMM_95

  !========================================
  ! HEMM implementation


  subroutine cHEMM_95(a,b,c,side,uplo,alpha,beta)
    complex(skind),dimension(:,:),intent(in):: a
    complex(skind),dimension(:,:),intent(in):: b
    complex(skind),dimension(:,:),intent(inout):: c
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    complex(skind),intent(in),optional:: beta
    complex(skind):: vbeta
    integer:: lda,ldb,ldc
    integer:: m,n,k

    vside = "L"
    if (present(side)) vside = side%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = c_one
    if (present(alpha)) valpha = alpha

    vbeta = c_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/csize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/csize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/csize
    m = size(c,1)
    n = size(c,2)
    call cHEMM(vside,vuplo,m,n,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine cHEMM_95

  subroutine zHEMM_95(a,b,c,side,uplo,alpha,beta)
    complex(dkind),dimension(:,:),intent(in):: a
    complex(dkind),dimension(:,:),intent(in):: b
    complex(dkind),dimension(:,:),intent(inout):: c
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    complex(dkind),intent(in),optional:: beta
    complex(dkind):: vbeta
    integer:: lda,ldb,ldc
    integer:: m,n,k

    vside = "L"
    if (present(side)) vside = side%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    valpha = z_one
    if (present(alpha)) valpha = alpha

    vbeta = z_zero
    if (present(beta)) vbeta = beta
    lda = 1
    if (size(a,2) > 1 .and. size(a,1) > 0) lda = (loc(a(1,2))-loc(a(1,1)))/zsize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/zsize
    ldc = 1
    if (size(c,2) > 1 .and. size(c,1) > 0) ldc = (loc(c(1,2))-loc(c(1,1)))/zsize
    m = size(c,1)
    n = size(c,2)
    call zHEMM(vside,vuplo,m,n,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine zHEMM_95

  !========================================
  ! TRMM implementation



  subroutine sTRMM_95(t,b,side,uplo,transt,diag,alpha)
    real(skind),dimension(:,:),intent(in):: t
    real(skind),dimension(:,:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    integer:: ldt,ldb
    integer:: m,n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c

    valpha = s_one
    if (present(alpha)) valpha = alpha
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/ssize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/ssize
    m = size(b,1)
    n = size(b,2)
    call sTRMM(vside,vuplo,vtranst,vdiag,m,n,valpha,t(1,1),ldt,b(1,1),ldb)
  end subroutine sTRMM_95

  subroutine dTRMM_95(t,b,side,uplo,transt,diag,alpha)
    real(dkind),dimension(:,:),intent(in):: t
    real(dkind),dimension(:,:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    integer:: ldt,ldb
    integer:: m,n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c

    valpha = d_one
    if (present(alpha)) valpha = alpha
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/dsize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/dsize
    m = size(b,1)
    n = size(b,2)
    call dTRMM(vside,vuplo,vtranst,vdiag,m,n,valpha,t(1,1),ldt,b(1,1),ldb)
  end subroutine dTRMM_95

  subroutine cTRMM_95(t,b,side,uplo,transt,diag,alpha)
    complex(skind),dimension(:,:),intent(in):: t
    complex(skind),dimension(:,:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    integer:: ldt,ldb
    integer:: m,n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c

    valpha = c_one
    if (present(alpha)) valpha = alpha
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/csize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/csize
    m = size(b,1)
    n = size(b,2)
    call cTRMM(vside,vuplo,vtranst,vdiag,m,n,valpha,t(1,1),ldt,b(1,1),ldb)
  end subroutine cTRMM_95

  subroutine zTRMM_95(t,b,side,uplo,transt,diag,alpha)
    complex(dkind),dimension(:,:),intent(in):: t
    complex(dkind),dimension(:,:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    integer:: ldt,ldb
    integer:: m,n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c

    valpha = z_one
    if (present(alpha)) valpha = alpha
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/zsize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/zsize
    m = size(b,1)
    n = size(b,2)
    call zTRMM(vside,vuplo,vtranst,vdiag,m,n,valpha,t(1,1),ldt,b(1,1),ldb)
  end subroutine zTRMM_95

  !========================================
  ! TRSM implementation



  subroutine sTRSM_95(t,b,side,uplo,transt,diag,alpha)
    real(skind),dimension(:,:),intent(in):: t
    real(skind),dimension(:,:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    integer:: ldt,ldb
    integer:: m,n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c

    valpha = s_one
    if (present(alpha)) valpha = alpha
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/ssize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/ssize
    m = size(b,1)
    n = size(b,2)
    call sTRSM(vside,vuplo,vtranst,vdiag,m,n,valpha,t(1,1),ldt,b(1,1),ldb)
  end subroutine sTRSM_95

  subroutine dTRSM_95(t,b,side,uplo,transt,diag,alpha)
    real(dkind),dimension(:,:),intent(in):: t
    real(dkind),dimension(:,:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    integer:: ldt,ldb
    integer:: m,n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c

    valpha = d_one
    if (present(alpha)) valpha = alpha
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/dsize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/dsize
    m = size(b,1)
    n = size(b,2)
    call dTRSM(vside,vuplo,vtranst,vdiag,m,n,valpha,t(1,1),ldt,b(1,1),ldb)
  end subroutine dTRSM_95

  subroutine cTRSM_95(t,b,side,uplo,transt,diag,alpha)
    complex(skind),dimension(:,:),intent(in):: t
    complex(skind),dimension(:,:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    integer:: ldt,ldb
    integer:: m,n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c

    valpha = c_one
    if (present(alpha)) valpha = alpha
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/csize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/csize
    m = size(b,1)
    n = size(b,2)
    call cTRSM(vside,vuplo,vtranst,vdiag,m,n,valpha,t(1,1),ldt,b(1,1),ldb)
  end subroutine cTRSM_95

  subroutine zTRSM_95(t,b,side,uplo,transt,diag,alpha)
    complex(dkind),dimension(:,:),intent(in):: t
    complex(dkind),dimension(:,:),intent(inout):: b
    type(blas_side_type),intent(in),optional:: side
    character:: vside
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: transt
    character:: vtranst
    type(blas_diag_type),intent(in),optional:: diag
    character:: vdiag
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    integer:: ldt,ldb
    integer:: m,n

    vside = "L"
    if (present(side)) vside = side%c

    vtranst = "N"
    if (present(transt)) vtranst = transt%c

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vdiag = "N"
    if (present(diag)) vdiag = diag%c

    valpha = z_one
    if (present(alpha)) valpha = alpha
    ldt = 1
    if (size(t,2) > 1 .and. size(t,1) > 0) ldt = (loc(t(1,2))-loc(t(1,1)))/zsize
    ldb = 1
    if (size(b,2) > 1 .and. size(b,1) > 0) ldb = (loc(b(1,2))-loc(b(1,1)))/zsize
    m = size(b,1)
    n = size(b,2)
    call zTRSM(vside,vuplo,vtranst,vdiag,m,n,valpha,t(1,1),ldt,b(1,1),ldb)
  end subroutine zTRSM_95

  !========================================
  ! SYRK implementation


  subroutine sSYRK_95(a,c,uplo,trans,alpha,beta)
    real(skind),dimension(:,:),intent(in):: a
    real(skind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    real(skind),intent(in),optional:: beta
    real(skind):: vbeta
    integer:: lda,ldc
    integer:: n,k

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = s_one
    if (present(alpha)) valpha = alpha

    vbeta = s_zero
    if (present(beta)) vbeta = beta
    n = size(c,1)
    if (vtrans == "N") then
       k = size(a,2)
    else
       k = size(a,1)
    end if
    call sSYRK(vuplo,vtrans,n,k,valpha,a(1,1),lda,vbeta,c(1,1),ldc)
  end subroutine sSYRK_95

  subroutine dSYRK_95(a,c,uplo,trans,alpha,beta)
    real(dkind),dimension(:,:),intent(in):: a
    real(dkind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    real(dkind),intent(in),optional:: beta
    real(dkind):: vbeta
    integer:: lda,ldc
    integer:: n,k

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = d_one
    if (present(alpha)) valpha = alpha

    vbeta = d_zero
    if (present(beta)) vbeta = beta
    n = size(c,1)
    if (vtrans == "N") then
       k = size(a,2)
    else
       k = size(a,1)
    end if
    call dSYRK(vuplo,vtrans,n,k,valpha,a(1,1),lda,vbeta,c(1,1),ldc)
  end subroutine dSYRK_95

  subroutine cSYRK_95(a,c,uplo,trans,alpha,beta)
    complex(skind),dimension(:,:),intent(in):: a
    complex(skind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    complex(skind),intent(in),optional:: beta
    complex(skind):: vbeta
    integer:: lda,ldc
    integer:: n,k

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = c_one
    if (present(alpha)) valpha = alpha

    vbeta = c_zero
    if (present(beta)) vbeta = beta
    n = size(c,1)
    if (vtrans == "N") then
       k = size(a,2)
    else
       k = size(a,1)
    end if
    call cSYRK(vuplo,vtrans,n,k,valpha,a(1,1),lda,vbeta,c(1,1),ldc)
  end subroutine cSYRK_95

  subroutine zSYRK_95(a,c,uplo,trans,alpha,beta)
    complex(dkind),dimension(:,:),intent(in):: a
    complex(dkind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    complex(dkind),intent(in),optional:: beta
    complex(dkind):: vbeta
    integer:: lda,ldc
    integer:: n,k

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = z_one
    if (present(alpha)) valpha = alpha

    vbeta = z_zero
    if (present(beta)) vbeta = beta
    n = size(c,1)
    if (vtrans == "N") then
       k = size(a,2)
    else
       k = size(a,1)
    end if
    call zSYRK(vuplo,vtrans,n,k,valpha,a(1,1),lda,vbeta,c(1,1),ldc)
  end subroutine zSYRK_95

  !========================================
  ! HERK implementation


  subroutine cHERK_95(a,c,uplo,trans,alpha,beta)
    complex(skind),dimension(:,:),intent(in):: a
    complex(skind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    complex(skind),intent(in),optional:: beta
    complex(skind):: vbeta
    integer:: lda,ldc
    integer:: n,k

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = c_one
    if (present(alpha)) valpha = alpha

    vbeta = c_zero
    if (present(beta)) vbeta = beta
    n = size(c,1)
    if (vtrans == "N") then
       k = size(a,2)
    else
       k = size(a,1)
    end if
    call cHERK(vuplo,vtrans,n,k,valpha,a(1,1),lda,vbeta,c(1,1),ldc)
  end subroutine cHERK_95

  subroutine zHERK_95(a,c,uplo,trans,alpha,beta)
    complex(dkind),dimension(:,:),intent(in):: a
    complex(dkind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    complex(dkind),intent(in),optional:: beta
    complex(dkind):: vbeta
    integer:: lda,ldc
    integer:: n,k

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = z_one
    if (present(alpha)) valpha = alpha

    vbeta = z_zero
    if (present(beta)) vbeta = beta
    n = size(c,1)
    if (vtrans == "N") then
       k = size(a,2)
    else
       k = size(a,1)
    end if
    call zHERK(vuplo,vtrans,n,k,valpha,a(1,1),lda,vbeta,c(1,1),ldc)
  end subroutine zHERK_95

  !========================================
  ! SYR2K implementation


  subroutine sSYR2K_95(a,b,c,uplo,trans,alpha,beta)
    real(skind),dimension(:,:),intent(in):: a,b
    real(skind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    real(skind),intent(in),optional:: alpha
    real(skind):: valpha
    real(skind),intent(in),optional:: beta
    real(skind):: vbeta
    integer:: lda,ldb,ldc
    integer:: n,k

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = s_one
    if (present(alpha)) valpha = alpha

    vbeta = s_zero
    if (present(beta)) vbeta = beta
    n = size(c,1)
    if (vtrans == "N") then
       k = size(a,2)
    else
       k = size(a,1)
    end if
    call sSYR2K(vuplo,vtrans,n,k,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine sSYR2K_95

  subroutine dSYR2K_95(a,b,c,uplo,trans,alpha,beta)
    real(dkind),dimension(:,:),intent(in):: a,b
    real(dkind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    real(dkind),intent(in),optional:: alpha
    real(dkind):: valpha
    real(dkind),intent(in),optional:: beta
    real(dkind):: vbeta
    integer:: lda,ldb,ldc
    integer:: n,k

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = d_one
    if (present(alpha)) valpha = alpha

    vbeta = d_zero
    if (present(beta)) vbeta = beta
    n = size(c,1)
    if (vtrans == "N") then
       k = size(a,2)
    else
       k = size(a,1)
    end if
    call dSYR2K(vuplo,vtrans,n,k,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine dSYR2K_95

  subroutine cSYR2K_95(a,b,c,uplo,trans,alpha,beta)
    complex(skind),dimension(:,:),intent(in):: a,b
    complex(skind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    complex(skind),intent(in),optional:: beta
    complex(skind):: vbeta
    integer:: lda,ldb,ldc
    integer:: n,k

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = c_one
    if (present(alpha)) valpha = alpha

    vbeta = c_zero
    if (present(beta)) vbeta = beta
    n = size(c,1)
    if (vtrans == "N") then
       k = size(a,2)
    else
       k = size(a,1)
    end if
    call cSYR2K(vuplo,vtrans,n,k,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine cSYR2K_95

  subroutine zSYR2K_95(a,b,c,uplo,trans,alpha,beta)
    complex(dkind),dimension(:,:),intent(in):: a,b
    complex(dkind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    complex(dkind),intent(in),optional:: beta
    complex(dkind):: vbeta
    integer:: lda,ldb,ldc
    integer:: n,k

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = z_one
    if (present(alpha)) valpha = alpha

    vbeta = z_zero
    if (present(beta)) vbeta = beta
    n = size(c,1)
    if (vtrans == "N") then
       k = size(a,2)
    else
       k = size(a,1)
    end if
    call zSYR2K(vuplo,vtrans,n,k,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine zSYR2K_95

  !========================================
  ! HER2K implementation


  subroutine cHER2K_95(a,b,c,uplo,trans,alpha,beta)
    complex(skind),dimension(:,:),intent(in):: a,b
    complex(skind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    complex(skind),intent(in),optional:: alpha
    complex(skind):: valpha
    complex(skind),intent(in),optional:: beta
    complex(skind):: vbeta
    integer:: lda,ldb,ldc
    integer:: n,k

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = c_one
    if (present(alpha)) valpha = alpha

    vbeta = c_zero
    if (present(beta)) vbeta = beta
    n = size(c,1)
    if (vtrans == "N") then
       k = size(a,2)
    else
       k = size(a,1)
    end if
    call cHER2K(vuplo,vtrans,n,k,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine cHER2K_95

  subroutine zHER2K_95(a,b,c,uplo,trans,alpha,beta)
    complex(dkind),dimension(:,:),intent(in):: a,b
    complex(dkind),dimension(:,:),intent(inout):: c
    type(blas_uplo_type),intent(in),optional:: uplo
    character:: vuplo
    type(blas_trans_type),intent(in),optional:: trans
    character:: vtrans
    complex(dkind),intent(in),optional:: alpha
    complex(dkind):: valpha
    complex(dkind),intent(in),optional:: beta
    complex(dkind):: vbeta
    integer:: lda,ldb,ldc
    integer:: n,k

    vuplo = "U"
    if (present(uplo)) vuplo = uplo%c

    vtrans = "N"
    if (present(trans)) vtrans = trans%c

    valpha = z_one
    if (present(alpha)) valpha = alpha

    vbeta = z_zero
    if (present(beta)) vbeta = beta
    n = size(c,1)
    if (vtrans == "N") then
       k = size(a,2)
    else
       k = size(a,1)
    end if
    call zHER2K(vuplo,vtrans,n,k,valpha,a(1,1),lda,b(1,1),ldb,vbeta,c(1,1),ldc)
  end subroutine zHER2K_95

end module blas95
