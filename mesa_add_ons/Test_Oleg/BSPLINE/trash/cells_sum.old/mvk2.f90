!=========================================================================
    SUBROUTINE mvk2(k)
!=========================================================================
!
!   vk matrix elements in the Spline basis with cell algorithm
!
!   Calling sequence:
!
!         mvk2
!       ---------
!      /        \\
!   momentd   vk_pdiag
!                  ||
!               vk_triang
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
!      rkb    four-dimensional array of vk integrals of power k in
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
    if(ntype == 'vk2' .and. knk == k) Return

    ! .. compute the moments in the spline basis

    CALL moments(  k+1 , rkd1)
    CALL moments(-(k+2), rkd2)
    CALL momentt(  k   , rkd3)
    CALL momentt(-(k+3), rkd4)
    CALL vk_pdiag

    ntype = 'vk2'
    knk = k

!-----------------------------------------------------------------------
    CONTAINS
!-----------------------------------------------------------------------

!======================================================================
    SUBROUTINE vk_pdiag
!======================================================================
!
!   Computes the vk matrix elements in the triangle cells
!
!   SUBROUTINES called:
!       vk_triang
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

    INTEGER :: iv
    REAL(KIND=8) :: hp1

    ! .. the first equal step region.

    Do iv=1,ml+ks-1
      CALL vk_triang(iv)
    End do

    ! .. the log region --- using scaling law.

    hp1=h+1.d0
    Do iv=ml+ks,ml+me-ks+2
      rkd(:,:,iv) = rkd(:,:,iv-1) / hp1
    End do

    ! .. the last equal step region

    DO iv=ml+me-ks+3,nv
      CALL vk_triang(iv)
    END DO


    END SUBROUTINE vk_pdiag



!========================================================================
    SUBROUTINE vk_triang(iv)
!========================================================================
!
!    Returns the "VK matrix element" in the diagonal cell   iv
!
!------------------------------------------------------------------------
!
!   SUBROUTINES called: gauss, bsplvd
!
!---------------------------------------------------------------------
!
!   On entry    iv  - index of the diagonal cell 
!   --------
!    
!       
!
!---------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: iv

! .. local variables

    INTEGER(4) :: i,j, ip,jp, m, left, ii,jj
    REAL(8) :: xbase, c
    REAL(8), DIMENSION(ks) :: x,w,gx,gw, bi
    REAL(8), DIMENSION(ks,ks) :: bspTmp,bspdTmp
    REAL(8), DIMENSION(ks,ks,ks) ::Int1,Int2
    REAL(8), DIMENSION(nv,ks,ks) :: dbiatx	          

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
! .. Int1  =  |      bsp(iv,:,i) bsp(iv,:,ip) r  dr
!             / r_iv

 
      c = grm(iv,m)**(k+2) * grw(iv,m)
      DO i=1,ks
       bi(:) = gw(:)*bspTmp(:,i)
       DO ip=1,ks
        Int1(i,ip,m)= SUM(bi(:)*bspdTmp(:,ip)) * c
       END DO
      END DO

!              / r(iv,m)                       k+1
! .. Int2  =  |      bsp(iv,:,j) bsp(iv,:,jp) r    dr
!             / r_iv

      gw = gw * gx
      c = grm(iv,m)**(k+3) * grw(iv,m)
      DO j=1,ks
       bi(:) = gw(:)*bspTmp(:,j)
       DO jp=1,ks
        Int2(j,jp,m)= SUM(bi(:)*bspTmp(:,jp)) * c
       END DO
      END DO

    END DO	!  over m

    jj = 0
    DO j=1,ks
     DO jp=1,ks
     jj = jj + 1
     gx(:) = bsp(iv,:,j)*bsp(iv,:,jp)

      ii = 0
      DO i=1,ks
       DO ip=1,ks
        ii = ii + 1

        rkd(ii,jj,iv) = SUM( gx(:)*Int1(i,ip,:))

       END DO
      END DO
     END DO
    END DO

      ii = 0
      DO i=1,ks
       DO ip=1,ks
        ii = ii + 1
        gx(:) = bsp(iv,:,i)*bsq(iv,:,ip)
   
        jj = 0
        DO j=1,ks
         DO jp=1,ks
          jj = jj + 1
          rkd(ii,jj,iv) = rkd(ii,jj,iv) + SUM( gx(:)*Int2(j,jp,:)) 

       END DO
      END DO
     END DO
    END DO

    END SUBROUTINE vk_triang

    END SUBROUTINE mvk2
