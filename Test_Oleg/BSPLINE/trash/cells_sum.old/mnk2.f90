!=========================================================================
    SUBROUTINE mnk2(k)
!=========================================================================
!
!            //           r2^k
!   Computes || dr1 dr2  -----B_i(r1) B_j(r2) B_ip(r1) B_jp(r2) E(r1-r2)
!            //         r1^(k+3)
!
!
!   Calling sequence:
!
!       mnk2
!       ----
!      /    \
!   moment nk_pdiag
!            ||
!          nk_triang
!           /   \
!        gauss  qbsplvb
!
!-------------------------------------------------------------------------
!
!   on entry
!   --------
!       k     the power of the factor in the integrals
!
!   on exit
!   -------
!       rkb   four-dimensional array of the Nk integrals in
!             the Spline basis
!
!-------------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_momentc
    USE spline_integrals

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k

    ! .. check the need of calculations

    if(ntype == 'aaa') Call allocate_momentc
    if(ntype == 'nk2' .and. knk == k) Return

    ! .. compute the moments in the spline basis
  
    CALL moments(  k   , rkd1)
    CALL moments(-(k+3), rkd2)
    CALL nk_pdiag
 
    ntype = 'nk2'
    knk = k

!-------------------------------------------------------------------
  CONTAINS
!-------------------------------------------------------------------

!===================================================================
  SUBROUTINE nk_pdiag
!===================================================================
!
!   Computes the Nk matrix elements in the triangle cells
!
!   SUBROUTINES called:
!       nk_triang
!
!-------------------------------------------------------------------
!
!   on entry
!   --------
!       k       the power of moments
!
!   on exit
!   -------
!       rkt     the four-dimensional array of pieces
!               defining <B_i B_j|r2^k/r1^(k+3) E(r1-r2)|B_i' B_j'>
!               over the lower triangle
!
!-------------------------------------------------------------------

    IMPLICIT NONE

    ! .. local variables

    INTEGER :: iv,ik
    REAL(KIND=8) :: hp1

    ! .. the first equal step region.

    DO iv=1,ml+ks-1
      CALL nk_triang(iv)
    END DO

    ! .. the log region --- using scaling law.

    hp1=h+1.d0
    ik = ks*(ks+1)/2
    DO iv=ml+ks,ml+me-ks+2
      rkd(1:ik,1:ik,iv) = rkd(1:ik,1:ik,iv-1) / hp1
    END DO

    ! .. the last equal step region

    DO iv=ml+me-ks+3,nv
      CALL nk_triang(iv)
    END DO

    END SUBROUTINE nk_pdiag


!========================================================================
    SUBROUTINE nk_triang(iv)
!========================================================================
!
!   Returns the "NK matrix element" in low triangle cell 
!                                                           
!------------------------------------------------------------------------
!
!   SUBROUTINES called:
!       gauss
!       vbsplvd
!
!------------------------------------------------------------------------
!
!   On entry
!   --------
!       k:     the indices of of the bsplines
!       iv:    the index of the integration region
!
!------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: iv

    ! .. local variables

    INTEGER :: i,j, ip,jp, m,left, ii,jj
    REAL(KIND=8) :: xbase, c
    REAL(KIND=8), DIMENSION(ks) :: x,w, gx,gw
    REAL(KIND=8), DIMENSION(ks,ks) :: bspTmp
    REAL(KIND=8), DIMENSION(ks,ks,ks) :: Int
    REAL(KIND=8), DIMENSION(nv,ks,ks) :: dbiatx

    left=iv+ks-1
    xbase=t(left)

    ! .. setup the gaussian points

    CALL gauss(ks,x,w)

    DO m=1,ks

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks
        Call vbsplvd(t,left,1,gx(i),1,dbiatx)
        bspTmp (i,:) = dbiatx(1,:,1)
      END DO

! .. and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:) * gx(:)**k

!            / r(iv,m)                             k
! .. Int =  |      bsp(iv,:,j)(r) bsp(iv,:,jp)(r) r  dr
!           / r_iv

      c = grw(iv,m) * grm(iv,m)**(k+3)
      DO j=1,ks
        gx(:) = gw(:)*bspTmp(:,j)
        DO jp=j,ks
          Int(j,jp,m) = SUM(gx(:)*bspTmp(:,jp)) * c
        END DO
      END DO

    END DO    ! over m

! .. second integration ..

    ii = 0
    DO i=1,ks
     DO ip=i,ks
      ii = ii + 1
      gx(:) = bsp(iv,:,i)*bsp(iv,:,ip)
      jj = 0
      DO j=1,ks
       DO jp=j,ks
         jj = jj + 1
         rkd(ii,jj,iv) = SUM(gx(:)*INT(j,jp,:))
       END DO
      END DO
     END DO
    END DO

    END SUBROUTINE nk_triang

    END SUBROUTINE mnk2

