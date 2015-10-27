!deck chk_nrm_08_1d.f
!***begin prologue     chk_nrm_08_1d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for one dimension wavefunction
!***
!***references
!***routines called    sdot
!***end prologue       chk_nrm_08_1d
  SUBROUTINE chk_nrm_08_1d(v)
  USE dvrprop_global
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: v
  REAL*8                                 :: ddot
  REAL*8                                 :: nrm, first, last
  INTEGER                                :: i
!
  nrm = ddot(nphy(1),v(1,1),1,v(1,1),1)   &
                        +                 &
        ddot(nphy(1),v(1,2),1,v(1,2),1)
  IF( typke == 'fd' ) then
      first = v(1,1) * v(1,1) + v(1,2) * v(1,2)
      last = v(nphy(1),1) * v(nphy(1),1) + v(nphy(1),2) * v(nphy(1),2)
      nrm = nrm - .5d0 * ( first + last)
      nrm = del * nrm
  END IF
  write(iout,1) nrm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE chk_nrm_08_1d
