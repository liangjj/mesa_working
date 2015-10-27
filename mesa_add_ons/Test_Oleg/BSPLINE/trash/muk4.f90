!=========================================================================
    SUBROUTINE muk4(k)
!=========================================================================
!
!   uk matrix elements in the Spline basis with cell algorithm
!
!   Calling sequence:
!
!         muk4
!       ---------
!      /        \\
!   momentd   uk_pdiag
!                  ||
!               uk_triang
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
!      rkd    the cell integrals for UK operator
!-------------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_galerkin
    USE spline_moments

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k

    ! .. check the need of calculations

    if(mtype == 'aaa') Call allocate_moments
    if(mtype == 'uk4' .and. kmk == k) Return

    ! .. compute the moments in the spline basis
    
    CALL moments(-(k+2), rkm1)
    CALL momentt(  k   , rkm2)
    CALL momentt(  k-1 , rkm3)
    CALL moments(-(k+1), rkm4)
    CALL uk_pdiag

    mtype = 'uk4'
    kmk = k

!-----------------------------------------------------------------------
    CONTAINS
!-----------------------------------------------------------------------

!======================================================================
    SUBROUTINE uk_pdiag
!======================================================================
!
!   Computes the uk matrix elements in the triangle cells
!
!   SUBROUTINES called:
!       uk_triang
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
!                  defining UK integral over diagonal cell
!
!----------------------------------------------------------------------

    IMPLICIT NONE

    ! .. local variables

    INTEGER :: iv, ik
    REAL(KIND=8) :: hp1

    ! .. the first equal step region.

    Do iv=1,ml+ks-1
      CALL uk_triang(iv)
    End do

    ! .. the log region --- using scaling law.

    ik = ks*(ks+1)/2
    hp1=h+1.d0
    Do iv=ml+ks,ml+me-ks+2
      rkd(1:ik,:,iv) = rkd(1:ik,:,iv-1) / hp1
    End do

    ! .. the last equal step region

    DO iv=ml+me-ks+3,nv
      CALL uk_triang(iv)
    END DO


    END SUBROUTINE uk_pdiag

!========================================================================
    SUBROUTINE uk_triang(iv)
!========================================================================
!
!    Returns the "uk matrix element" in the diagonal cell "iv"
!
!------------------------------------------------------------------------
!
!   SUBROUTINES called:
!       gauss
!       bsplvd
!
!------------------------------------------------------------------------
!
!   On entry
!   --------
!       k:     the indices of of the bsplines
!       iv:    the index of the integration region
!
!   On exit
!   -------
!      rkd:    array of diagonal cell integrals
!------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: iv

! .. local variables

    INTEGER :: i, j, ip, jp, m, left, ik, ii, jj
    REAL(KIND=8) :: xbase, c
    REAL(KIND=8), DIMENSION(ks) :: x,w, gx,gw, bi
    REAL(KIND=8), DIMENSION(ks,ks) :: bspTmp,bspdTmp
    REAL(KIND=8), DIMENSION(ks,ks*ks) ::Int1,Int2
    REAL(KIND=8), DIMENSION(ks*ks,ks*ks) ::a,b
    REAL(KIND=8), DIMENSION(nv,ks,ks) :: dbiatx	

    ik = ks*(ks+1)/2
    left = iv+ks-1
    xbase = t(left)

! .. setup the gaussian points

    CALL gauss(ks,x,w)

    DO m=1,ks                  !  loop over old knots

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks               
       Call vbsplvd(t,left,1,gx(i),2,dbiatx)
       bspTmp (i,1:ks) = dbiatx(1,1:ks,1)
       bspdTmp(i,1:ks) = dbiatx(1,1:ks,2) - bspTmp(i,1:ks)/gx(i)
      END DO

! .. and the corresponding gaussian weights
      
      gw(:) = (gr(iv,m)-xbase)*w(:) * gx(:)**k


!              / r(iv,m)         ___           k
! .. Int1  =  |      bsp(iv,:,j) bsp(iv,:,jp) r  dr
!             / r_iv

 
      c = grm(iv,m)**(k+2) * grw(iv,m) * (k-1)
      jj = 0
      DO j=1,ks
       bi(:) = gw(:)*bspTmp(:,j)
       DO jp=1,ks
        jj = jj + 1
        Int1(m,jj)= SUM(bi(:)*bspdTmp(:,jp)) * c
       END DO
      END DO

!              / r(iv,m)                       k-1
! .. Int2  =  |      bsp(iv,:,i) bsp(iv,:,ip) r    dr
!             / r_iv

      gw = gw / gx
      c = grm(iv,m)**(k+1) * grw(iv,m) * (k+2)
      ii = 0
      DO i=1,ks
       bi(:) = gw(:)*bspTmp(:,i)
       DO ip=i,ks
        ii = ii + 1
        Int2(m,ii)= SUM(bi(:)*bspTmp(:,ip)) * c
       END DO
      END DO

    END DO	!  over m

! ..  second integration ..

      Do i=1,ks
       bspTmp (:,i) = bsp(iv,:,i)
       bspdTmp(:,i) = bspd(iv,:,i,1) - grm(iv,:)*bsp(iv,:,i)
      End do

      ii = 0
      DO i=1,ks
       DO ip=i,ks
        ii = ii + 1
        bi(:) = bspTmp(:,i)*bspTmp(:,ip)
        jj = 0
        DO j=1,ks
         DO jp=1,ks
          jj = jj + 1
          b(ii,jj) = SUM(bi(:)*Int1(:,jj)) 
         END DO
        END DO
       END DO
      END DO

      jj = 0
      DO j=1,ks
       DO jp=1,ks
        jj = jj + 1
        bi (:) = bspTmp(:,j)*bspdTmp(:,jp)
        ii = 0
        DO i=1,ks
         DO ip=i,ks
          ii = ii + 1
          a(ii,jj) = SUM(bi(:)*Int2(:,ii))
         END DO
        END DO
       END DO
      END DO
     
    rkd(1:ik,:,iv) =  a(1:ik,:) + b(1:ik,:)

    END SUBROUTINE uk_triang

    END SUBROUTINE muk4
