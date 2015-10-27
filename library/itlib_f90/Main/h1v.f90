!deck h1v.f
!***begin prologue     h1v
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            space hamiltonian times vector
!***
!***references
 
!***routines called
!***end prologue       h1v
!
  SUBROUTINE h1v(vin,vout,nvc)
  USE dvd_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                                :: nvc
  REAL*8, DIMENSION(n_dvd,nvc)           :: vin, vout
  CALL apbc(vout,grid(1)%h,vin,nphy(1),nphy(1),nvc)
END SUBROUTINE h1v

