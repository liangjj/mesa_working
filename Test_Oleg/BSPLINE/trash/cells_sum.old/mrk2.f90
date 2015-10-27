!=========================================================================
    SUBROUTINE mrk2(k)
!=========================================================================
!
!   Defines moments for Rk-integrals in the B-spline cells
!
!   Calling sequence:            mrk2             
!                                ----             
!                               /    \\           
!                           moments rk_pdiag      
!                                     ||          
!                                   rk_triang     
!                                    /   \        
!                                 gauss  vbsplvb  
! 
!-------------------------------------------------------------------------
!
!   on entry    k  -  multipole index
!   --------
!       
!   on exit     rkd1,rkd2,rkd - off-diagonal and diagonal moments 
!   -------                     in the symmetric mode
!
!-------------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_momentc

    IMPLICIT NONE
    INTEGER(4), INTENT(IN) :: k

    ! .. check the need of calculations

    if(ntype == 'aaa') Call allocate_momentc
    if(ntype == 'rk2' .and. knk == k) Return

    ! .. compute the moments in the spline basis

    CALL moments(  k   , rkd1)
    CALL moments(-(k+1), rkd2)
    CALL rk_pdiag

    ntype='rk2'
    knk=k

    CONTAINS


!======================================================================
    SUBROUTINE rk_pdiag
!======================================================================
!
!   Computes the Slater matrix elements in the diagonal cells
!
!   SUBROUTINES called:  rk_triang
!
!----------------------------------------------------------------------
!
!   on entry       k  - multipole index
!   --------
!       
!   on exit        rkd - the four-dimensional array of pieces
!   -------              defining <B_i B_j|r^k/r^(k+1)|B_i' B_j'>
!                        over a triangle in given diagonal cell
!                  
!----------------------------------------------------------------------

    IMPLICIT NONE

    ! .. local variables

    INTEGER(4) :: iv,ik
    REAL(8) :: hp1

    ! .. the first equal step region.

    DO iv=1,ml+ks-1
      CALL rk_triang(iv)
    END DO

    ! .. the log region --- using scaling law.

    hp1=h+1.d0
    ik = ks*(ks+1)/2
    DO iv=ml+ks,ml+me-ks+2
      rkd(1:ik,1:ik,iv) =  rkd(1:ik,1:ik,iv-1) * hp1
    END DO

    ! .. the last equal step region

    DO iv=ml+me-ks+3,nv
      CALL rk_triang(iv)
    END DO

    END SUBROUTINE rk_pdiag



!========================================================================
    SUBROUTINE rk_triang(iv)
!========================================================================
!
!                                                            iv+1 -- iv+1
!   Returns the "slater matrix element" in triangle cell       |       |  ,
!                                                            iv   -- iv+1
!   i. e.,
!
!  /
! /              k
! |             r_<
! | dr_1 dr_2   ---  bsp(iv,:,i)(r_1) bsp(iv,:,j)(r_2)
! |              k+1
! |             r_>
! /                  bsp(iv,:,ip)(r_1) bsp(iv,:,jp)(r_2)
!/
!------------------------------------------------------------------------
!
!   SUBROUTINES called:   gauss, vbsplvd
!
!---------------------------------------------------------------------
!
!   On entry   iv  -  index of the interval (diag. cell)
!   --------
!
!   On exit    rkd(.,.,iv) - arrays of B-spline integrals for given 
!   --------                 interval iv in symmetric mode
!
!---------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: iv

    ! .. local variables

    INTEGER :: i, j, ip, jp, m, left, ii,jj, ik
    REAL(KIND=8) :: xbase
    REAL(KIND=8), DIMENSION(ks) :: x,w, bi, gx,gw
    REAL(KIND=8), DIMENSION(ks,ks) :: bspTmp
    REAL(KIND=8), DIMENSION(ks,ks,ks) :: INT
    REAL(KIND=8), DIMENSION(nv,ks,ks) :: dbiatx
    REAL(KIND=8), DIMENSION(ks*(ks+1)/2,ks*(ks+1)/2) :: a

    left=iv+ks-1
    xbase=t(left)

! .. setup the gaussian points

    CALL gauss(ks,x,w)

    DO m=1,ks

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks
       CALL vbsplvd(t,left,1,gx(i),1,dbiatx)
       bspTmp(i,1:ks)= dbiatx(1,1:ks,1)
      END DO

! .. and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:)

      IF(k>1) THEN
        gx(:) = gw(:)*gx(:)**k
      ELSE IF(k==1) THEN
        gx(:) = gw(:)*gx(:)
      ELSE IF(k==0) THEN
        gx(:) = gw(:)
      END IF

!            / r(iv,m)                             k
! .. INT =  |      bsp(iv,:,j)(r) bsp(iv,:,jp)(r) r  dr
!           / r_iv

      DO j=1,ks
       gw(:) = gx(:)*bspTmp(:,j)
       DO jp=j,ks
        INT(j,jp,m)= SUM(gw(:)*bspTmp(:,jp))
       END DO
      END DO
    
    END DO	 !  over m


    IF(k/=0) THEN
      gx(:) = grw(iv,:)*grm(iv,:)**(k+1)
    ELSE
      gx(:) = grw(iv,:)*grm(iv,:)
    END IF

    ii = 0
    DO i=1,ks
     DO ip=i,ks
     ii = ii+1

      bi(:) = bsp(iv,:,i)*bsp(iv,:,ip)*gx(:)

      jj = 0
      DO j=1,ks
       DO jp=j,ks
       jj = jj + 1

         a(ii,jj) =  SUM(bi(:)*INT(j,jp,:))

       END DO
      END DO
     END DO
    END DO
    
    ik = ks*(ks+1)/2
    rkd(1:ik,1:ik,iv) = a + TRANSPOSE(a)

    END SUBROUTINE rk_triang

    END SUBROUTINE mrk2
