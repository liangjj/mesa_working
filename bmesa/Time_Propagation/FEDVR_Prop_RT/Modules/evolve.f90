!deck evolve.f
!***begin prologue     evolve
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       evolve
  SUBROUTINE evolve(vector,vector_0,tau)
  USE arnoldi_global
  USE dvr_shared  
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: vector, vector_0
  REAL*8                                 :: tau
!
  call iosys('read real "initial state" from bec',         &
              2*nphy(1),vector_0,0,' ')
  vector_0(:,2) = - sin(energy*tau) * vector_0(:,1)
  vector_0(:,1) = cos(energy*tau) * vector_0(:,1)

  title='exact vector'
  call prntfmn(title,vector_0,nphy(1),2,nphy(1),2,iout,'e')
  title='computed vector'
  call prntfmn(title,vector,nphy(1),2,nphy(1),2,iout,'e')
END SUBROUTINE evolve
