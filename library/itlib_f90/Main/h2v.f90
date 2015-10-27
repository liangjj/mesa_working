!*deck h2v.f
!***begin prologue     h2v
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source             
!***purpose            space hamiltonian times space*time vector
!***                   
!***references         
!
!***routines called    
!***end prologue       h2v
  subroutine h2v(vin,vout,nvc)
  USE dvd_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                                :: nvc, i
  REAL*8, DIMENSION(nphy(2),nphy(1),nvc) :: vin, vout
  call apbc(vout,grid(2)%h,vin,nphy(2),nphy(2),nphy(1)*nvc)
  DO i=1,nvc
     call apbct(vout(1,1,i),vin(1,1,i),grid(1)%h,nphy(2),nphy(1),nphy(1))
  END DO
END SUBROUTINE h2v
