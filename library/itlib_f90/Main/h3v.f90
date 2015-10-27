!deck h3v.f
!***begin prologue     h3v
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            space hamiltonian times space*time vector
!***
!***references

!***routines called
!***end prologue       h3v
  SUBROUTINE h3v(vin,vout,nvc)
  USE dvd_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                                         :: nvc, i, j 
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nvc)  :: vin, vout     
  CALL apbc(vout,grid(3)%h,vin,nphy(3),nphy(3),nphy(2)*nphy(1)*nvc)
  DO  i=1,nphy(1)
      DO  j=1,nvc
          CALL apbct(vout(1,1,i,j),vin(1,1,i,j),grid(2)%h,nphy(3),  &
                                                       nphy(2),nphy(2))
      END DO
  END DO
  DO  i=1,nvc
      CALL apbct(vout(1,1,1,i),vin(1,1,1,i),grid(1)%h,              &
                                         nphy(3)*nphy(2),nphy(1),nphy(1))
  END DO
END SUBROUTINE h3v
