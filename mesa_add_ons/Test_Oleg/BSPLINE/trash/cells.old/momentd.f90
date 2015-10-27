!=====================================================================
    SUBROUTINE momentd (k,rkm)
!=====================================================================
!
!   Computes moments with derivatives in intervals
!
!---------------------------------------------------------------------
!
!   on entry
!   --------
!       k      the power of moments
!
!   on exit
!   -------
!                                               _
!       rkm    the moments defining as <B_i|r^k|B_j> over an interval
!
!---------------------------------------------------------------------

    USE spline_param
    USE spline_grid

    IMPLICIT NONE
    INTEGER, INTENT(in) :: k
    REAL(KIND=8), DIMENSION(ks,ks,nv), INTENT(out) :: rkm

    ! .. local variables

    INTEGER :: iv, i, j
    REAL(KIND=8) :: hp1
    REAL(KIND=8), DIMENSION(nv,ks) :: gw

    gw = grw
    if( k > 0 ) gw = gw * gr**k
    if( k < 0 ) gw = gw * grm**(-k)

      ! .. the first equal-step region

      DO iv=1,ml+ks-1
        DO i=1,ks
          DO j=1,ks
            rkm(i,j,iv) = SUM(bsp(iv,:,i)*gw(iv,:)*    &
                             (bspd(iv,:,j,1)-grm(iv,:)*bsp(iv,:,j)) )
          END DO
        END DO
      END DO

      ! .. the log region --- using scaling law

      hp1 = (1.d0+h)**k
      DO iv=ml+ks,ml+me-ks+2
          rkm(:,:,iv) = rkm(:,:,iv-1)*hp1
      END DO

      ! .. the last equal step region

      DO iv=ml+me-ks+3,nv
        DO i=1,ks
          DO j=1,ks
            rkm(i,j,iv) = SUM(bsp(iv,:,i)*gw(iv,:)*    &
                             (bspd(iv,:,j,1)-grm(iv,:)*bsp(iv,:,j)) )
          END DO
        END DO
      END DO

  END SUBROUTINE momentd

