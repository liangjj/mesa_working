!deck chk_nrm_08_2d.f
!***begin prologue     chk_nrm_08_2d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for two dimension wavefunction
!***
!***references
!***routines called    sdot
!***end prologue       chk_nrm_08_2d
  SUBROUTINE chk_nrm_08_2d(v,fac_2d)
  USE dvrprop_global
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: v
  REAL*8, DIMENSION(nphy(1))             :: fac_2d
  REAL*8                                 :: ddot
  REAL*8                                 :: nrm, first, last
  INTEGER                                :: i, j
!
  nrm = 0.d0
  DO i=1,nphy(1)
     fac_2d(i) = ddot(nphy(2),v(1,i,1),1,v(1,i,1),1)  &
                               +                   &
              ddot(nphy(2),v(1,i,2),1,v(1,i,2),1)
  END DO
  IF( typke == 'fd' ) then
      DO i=1,nphy(1)
         first = v(1,i,1) * v(1,i,1)                   &
                          +                            &
                 v(1,i,2) * v(1,i,2)
         last  = v(nphy(2),i,1) * v(nphy(2),i,1)       &
                                +                      &
                 v(nphy(2),i,2) * v(nphy(2),i,2)
         fac_2d(i) = fac_2d(i) - .5d0 * ( first + last)
      END DO
      fac_2d(:) = del * fac_2d(:)
  END IF
  DO i=1,nphy(1)
     nrm = nrm + fac_2d(i)
  END DO
  IF( typke == 'fd' ) then
      nrm = nrm - .5d0 * ( fac_2d(1) + fac_2d(nphy(1)) )
      nrm = del * nrm
  END IF
  write(iout,1) nrm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE chk_nrm_08_2d
