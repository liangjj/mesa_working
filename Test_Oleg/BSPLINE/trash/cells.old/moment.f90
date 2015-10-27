!=====================================================================
    SUBROUTINE moment (k,rkm)
!=====================================================================
!
!   Computes moments in intervals
!
!---------------------------------------------------------------------
!
!   on entry
!   --------
!       k       the power of moments
!
!   on exit
!   -------
!       rkm     the moments defining as <B_i|r^k|B_j> over an interval
!               (only for j > i )
!---------------------------------------------------------------------

    USE spline_param
    USE spline_grid

    IMPLICIT NONE

    INTEGER, INTENT(in) :: k
    REAL(8), DIMENSION(ks,ks,nv), INTENT(out) :: rkm

    ! .. local variables

    INTEGER(4) :: iv, i, j
    REAL(8) :: hp1, c
    REAL(8), DIMENSION(nv,ks) :: gw
    REAL(8), DIMENSION(ks) :: bi
    
      gw = grw
      if( k > 0 ) gw = gw * gr**k
      if( k < 0 ) gw = gw * grm**(-k)

      ! .. the initial equal-step region

      DO iv=1,ml+ks-1
        DO i=1,ks
          bi(:) = bsp(iv,:,i)*gw(iv,:)
          DO j=i,ks
            c = SUM(bi(:)*bsp(iv,:,j))
            rkm(i,j,iv) = c
            rkm(j,i,iv) = c
          END DO
        END DO
      END DO

      ! .. the log region --- using scaling law

      hp1= (h+1.d0)**(1+k)
      DO iv=ml+ks,ml+me-ks+2
        DO i=1,ks
          DO j=1,ks
            rkm(i,j,iv) = rkm(i,j,iv-1)*hp1
          END DO
        END DO
      END DO

      ! .. the last equal-step region

      DO iv=ml+me-ks+3,nv
        DO i=1,ks
          bi(:) = bsp(iv,:,i)*gw(iv,:)
          DO j=i,ks
            c = SUM(bi(:)*bsp(iv,:,j))
            rkm(i,j,iv) = c
            rkm(j,i,iv) = c
          END DO
        END DO
      END DO

  END SUBROUTINE moment

