!deck check_norm_1d.f
!***begin prologue     check_norm_1d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for one dimension wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_norm_1d
  SUBROUTINE check_norm_1d(v,norm)
  USE dvrprop_global_rt
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: v
  REAL*8                                 :: ddot
  REAL*8                                 :: norm, first, last
!
  norm = ddot(nphy(1),v(1,1),1,v(1,1),1)               &
                         +                             &
         ddot(nphy(1),v(1,2),1,v(1,2),1)
  IF( typke == 'fd' ) then
      first = v(1,1) * v(1,1) + v(1,2) * v(1,2)
      last  = v(nphy(1),1) * v(nphy(1),1)              &
                          +                            &
              v(nphy(1),2) * v(nphy(1),2)
      norm = norm - .5d0 * ( first + last)
      norm = del * norm
  END IF
  write(iout,1) norm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE check_norm_1d
