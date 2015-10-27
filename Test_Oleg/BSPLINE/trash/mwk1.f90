!=========================================================================
    SUBROUTINE mwk1(k)
!=========================================================================
!
!   Assembles wk matrix elements of power k in the Spline basis
!   with cell algorithm.
!
!   Calling sequence:
!
!              mwk1
!       -------------------
!      /        |         \\
!   moment  momentd    wk_pdiag
!                          ||
!                      wk_triang
!                        /      \
!                      gauss  vbsplvd
!
!-------------------------------------------------------------------------
!
!   on entry
!   --------
!       k     the multipole of the slater integrals
!
!   on exit
!   -------
!      rkb    four-dimensional array of wk integrals of power k in
!             the Spline basis
!
!-------------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_galerkin
    USE spline_integrals
    USE spline_moments
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k

    ! .. local variables

    INTEGER :: iv,ih,i, ihp,ip, jv,jh,j, jhp,jp
    REAL(8) :: c

    ! .. check the need of calculations

    if(itype == 'aaa') Call allocate_integrals
    if(itype == 'wk1' .and. krk == k) Return

    ! .. compute the moments in the spline basis

    if(mtype == 'aaa') Call allocate_moments
    
    CALL momentd(k  , rkmB)     
    CALL moment (k-1, rkm)
    rkmB = 2*rkmB + (k+2)*rkm
    CALL moment (-(k+2), rkm)  

    CALL momentd(-(k+3), rkkmB)  
    CALL moment (-(k+4), rkkm) 
    rkkmB = 2*rkkmB - (k+1)*rkkm    
    CALL moment (k+1, rkkm)      

    CALL wk_pdiag

    rkb = 0.d0

    DO jv=1,nv

      DO jh=1,ks
        j = jv + jh - 1
        DO jhp=jh,ks
          jp = jhp - jh + 1

          DO iv=1,nv

            DO ih=1,ks
              i = iv + ih -1
              DO ihp=1,ks
                ip = ihp - ih + ks

         IF( iv < jv ) THEN
                 
                  c =  rkmB(ih,ihp,iv)*rkm(jh,jhp,jv)
        
         ELSE IF( iv > jv ) THEN
                  
                  c =  rkkm(jh,jhp,jv)*rkkmB(ih,ihp,iv)
         ELSE
                  c =  rkt(ih,jh,ihp,jhp,iv)
         END IF

               rkb(i,j,ip,jp) = rkb(i,j,ip,jp) +  c * fine
              
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    itype = 'wk1'
    krk = k

!-----------------------------------------------------------------------
    CONTAINS
!-----------------------------------------------------------------------


!======================================================================
    SUBROUTINE wk_pdiag
!======================================================================
!            !   Computes the V^k matrix elements in the triangle cells
!
!   SUBROUTINES called:
!       wk_triang
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
!                                                _
!                  defining <B_i B_j|r^k/r^(k+3)|B_i' B_j'>
!                  over the lower triangle
!                  plus     <B_i B_j|r^k/r^(k+3)|B_i' B_j'>
!                  over the upper triangle
!
!----------------------------------------------------------------------

    IMPLICIT NONE

    ! .. local variables

    INTEGER :: iv
    REAL(KIND=8) :: hp1

    ! .. the first equal step region.

    DO iv=1,ml+ks-1
      CALL wk_triang(iv)
    END DO

    ! .. the log region --- using scaling law.

    hp1=h+1.d0
    DO iv=ml+ks,ml+me-ks+2
         rkt(:,:,:,:,iv) = rkt(:,:,:,:,iv-1) / hp1
    END DO

    ! .. the last equal step region

    DO iv=ml+me-ks+3,nv
      CALL wk_triang(iv)
    END DO

    END SUBROUTINE wk_pdiag


!========================================================================
    SUBROUTINE wk_triang(iv)
!========================================================================
!
!    Returns the "wk matrix element" in the diagonal cell   iv
!
!------------------------------------------------------------------------
!
!   SUBROUTINES called:
!       gauss
!       bsplvd
!
!---------------------------------------------------------------------
!
!   On entry
!   --------
!       k:     the indices of of the bsplines
!       iv:    the index of the integration region
!
!---------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(in) :: iv

! .. local variables

    INTEGER :: i, j, ip, jp, m, left
    REAL(KIND=8) :: xbase, c
    REAL(KIND=8), DIMENSION(ks) :: x,w,gx,gw,gv, bi
    REAL(KIND=8), DIMENSION(ks,ks) :: bspTmp,bspdTmp
    REAL(KIND=8), DIMENSION(ks,ks,ks) ::Int1,Int2
    REAL(KIND=8), DIMENSION(nv,ks,ks) :: dbiatx	          
    REAL(KIND=8), DIMENSION(ks,ks,ks,ks) ::a,b


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
       bspdTmp(i,1:ks) = 2*dbiatx(1,1:ks,2) + k*bspTmp(i,1:ks)/gx(i)
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

      gv = gw * gx
      c = grm(iv,m)**(k+3) * grw(iv,m)
      DO j=1,ks
       bi(:) = gv(:)*bspTmp(:,j)
       DO jp=j,ks
        Int2(j,jp,m)= SUM(bi(:)*bspTmp(:,jp)) * c
       END DO
      END DO


    END DO	!  over m


      DO j=1,ks
       DO jp=j,ks
        gx (:) = bsp(iv,:,j)*bsp(iv,:,jp)
        DO i=1,ks
         DO ip=1,ks
          a(i,j,ip,jp) = SUM( gx(:)*Int1(i,ip,:))
         END DO
        END DO
       END DO
      END DO

      Do ip=1,ks
       bspTmp(:,ip) = 2*bspd(iv,:,ip,1) - (k+3)*grm(iv,:)*bsp(iv,:,ip)
      End do

      DO i=1,ks
       DO ip=1,ks
        gx(:) = bsp(iv,:,i)*bspTmp(:,ip)
        DO j=1,ks
         DO jp=j,ks
          b(i,j,ip,jp) = SUM( gx(:)*Int2(j,jp,:)) 
         END DO
        END DO
       END DO
      END DO

     
      rkt(:,:,:,:,iv) =  a + b
   
   
    END SUBROUTINE wk_triang


    END SUBROUTINE mwk1
