!=========================================================================
    SUBROUTINE mtk2(k)
!=========================================================================
!
!   Tk matrix elements in the Spline basis with cell algorithm
!
!   Calling sequence:
!
!         mtk2
!       ---------
!      /        \\
!   momentd   tk_pdiag
!                  ||
!               tk_triang
!                /      \
!             gauss   vbsplvd
!
!-------------------------------------------------------------------------
!
!   on entry
!   --------
!       k     the power of the slater integrals
!
!   on exit
!   -------
!      rkd    four-dimensional array of Tk integrals of power k in
!             the Spline basis 
!-------------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_galerkin
    USE spline_momentc

    IMPLICIT NONE
    INTEGER, INTENT(in) :: k  

    ! .. check the need of calculations

    if(ntype == 'aaa') Call allocate_moments
    if(ntype == 'tk2' .and. knk == k) Return

    ! .. compute the moments in the spline basis

    CALL momentt(  k   , rkd1)
    CALL momentt(-(k+1), rkd2)

    CALL tk_pdiag

    ntype = 'tk2'
    knk = k

!-----------------------------------------------------------------------
    CONTAINS
!-----------------------------------------------------------------------

!======================================================================
    SUBROUTINE tk_pdiag
!======================================================================
!
!   Computes the Tk matrix elements in the triangle cells
!
!   SUBROUTINES called:
!       tk_triang
!
!----------------------------------------------------------------------
!
!   on entry
!   --------
!       k          the power of moments
!
!   on exit
!   -------
!       rkt        the four-dimensional array of pieces
!                                                _    _
!                  defining <B_i B_j|r^k/r^(k+1)|B_i' B_j'>
!
!                  over a triangle
!
!----------------------------------------------------------------------

    IMPLICIT NONE

    ! .. local variables

    INTEGER :: iv, jk
    REAL(KIND=8) :: hp1

    ! .. the first equal step region.

    Do iv=1,ml+ks-1
      CALL tk_triang(iv)
    End do

    ! .. the log region --- using scaling law.

    hp1=h+1.d0
    jk = ks*ks
    Do iv=ml+ks,ml+me-ks+2
      rkd(1:jk,1:jk,iv) = rkd(1:jk,1:jk,iv-1) / hp1
    End do

    ! .. the last equal step region

    DO iv=ml+me-ks+3,nv
      CALL tk_triang(iv)
    END DO


    END SUBROUTINE tk_pdiag



!========================================================================
    SUBROUTINE tk_triang(iv)
!========================================================================
!
!   Returns the "TK matrix element" in the diagonal cell  'iv'      
!                                                     
!
!  /
! /              k
! |             r_<
! | dr_1 dr_2   ---  bsp(iv,:,i)(r_1) bsp(iv,:,j)(r_2)
! |              k+1
! |             r_>  ___               ___
! /                  bsp(iv,:,ip)(r_1) bsp(iv,:,jp)(r_2)
!/
!------------------------------------------------------------------------
!
!   SUBROUTINES called:
!       gauss
!       vbsplvd
!
!---------------------------------------------------------------------
!
!   On entry
!   --------
!       k:     multipole index
!       iv:    the index of the integration region
!
!---------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(in) :: iv

! .. local variables

      INTEGER :: i,j, ii,jj, ip,jp, m, left, jk
      REAL(KIND=8) :: xbase
      REAL(KIND=8), DIMENSION(ks) :: x,w, gx,gw
      REAL(KIND=8), DIMENSION(ks,ks) :: bspTmp, bspdTmp
      REAL(KIND=8), DIMENSION(ks,ks,ks) :: INT
      REAL(KIND=8), DIMENSION(nv,ks,ks) :: dbiatx
      REAL(KIND=8), DIMENSION(ks*ks,ks*ks) :: a

      left=iv+ks-1
      xbase=t(left)

! .. setup the gaussian points

      CALL gauss(ks,x,w)

! .. first integration 

      DO m=1,ks

! .. the absolute coordinate at the new gaussian point

       gx(1:ks) = (gr(iv,m)-xbase)*x(1:ks) + xbase

! .. the bspline values at the new gaussian points

       DO i=1,ks
        Call vbsplvd(t,left,1,gx(i),2,dbiatx)
        bspTmp (i,1:ks) = dbiatx(1,1:ks,1)
        bspdTmp(i,1:ks) = dbiatx(1,1:ks,2)-bspTmp(i,1:ks)/gx(i)
      END DO

! .. and the corresponding gaussian weights

       gw(1:ks) = (gr(iv,m)-xbase)*w(1:ks)
       gx(1:ks) = gx(1:ks)**k * gw(1:ks)

!            / r(iv,m)         ___           k
! .. INT =  |      bsp(iv,:,j) bsp(iv,:,jp) r  dr
!           / r_iv

       DO j=1,ks
        gw(1:ks) = gx(1:ks) * bspTmp(1:ks,j)
        DO jp=1,ks
          INT(j,jp,m)= SUM(gw(1:ks) * bspdTmp(1:ks,jp))
        END DO
       END DO

      END DO    ! over m

! .. second integration

      gx(1:ks) = grw(iv,1:ks)*grm(iv,1:ks)**(k+1)

      Do ip=1,ks
       bspTmp(1:ks,ip) = bsq(iv,1:ks,ip)*gx(1:ks)
      End do

      ii = 0
      DO i=1,ks
       DO ip=1,ks
        ii = ii + 1
        gx(1:ks) =  bsp(iv,1:ks,i)*bspTmp(1:ks,ip)
        
        jj = 0
        DO j=1,ks
         DO jp=1,ks
          jj = jj + 1
          
          a(ii,jj) = SUM(gx(1:ks)*INT(j,jp,1:ks))

         END DO
        END DO
 
       END DO
      END DO

      jk = ks*ks
      rkd(1:jk,1:jk,iv) = a + Transpose(a)

      END SUBROUTINE tk_triang

      END SUBROUTINE mtk2
