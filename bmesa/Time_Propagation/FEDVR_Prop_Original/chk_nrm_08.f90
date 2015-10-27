!deck chk_nrm_08.f
!***begin prologue     chk_nrm_08
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization 
!***
!***references
!***routines called    sdot
!***end prologue       chk_nrm_08
  SUBROUTINE chk_nrm_08(v,n)
  USE dvrprop_global
  USE fd_global,         ONLY : del
  USE dvr_shared,        ONLY : typke
  USE dvr_global,        ONLY : spdim
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n,2)                 :: v
  REAL*8                                 :: sdot
  REAL*8                                 :: nrm
!
  IF( typke == 'fd' ) then
      nrm = sdot(n-2,v(2,1),1,v(2,1),1) + sdot(n-2,v(2,2),1,v(2,2),1)  &
            + .5d0 * ( v(1,1) * v(1,1) + v(n,1) * v(n,1) )             &
            + .5d0 * ( v(1,2) * v(1,2) + v(n,2) * v(n,2) )
      IF(spdim==1) THEN
         nrm = del * nrm
      ELSE IF(spdim==2) THEN
         nrm = del * del * nrm
      ELSE IF(spdim==3) THEN
         nrm = del * del * del * nrm
      END IF
  ELSE
      nrm = sdot(n3d,v(1,1),1,v(1,1),1) + sdot(n3d,v(1,2),1,v(1,2),1)
  END IF
  write(iout,1) nrm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE chk_nrm_08
