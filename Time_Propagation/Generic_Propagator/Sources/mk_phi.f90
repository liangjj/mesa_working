!deck mk_phi.f
!***begin prologue     mk_phi
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose
!***references
!***routines called
!***end prologue       mk_phi
  SUBROUTINE mk_phi(phi,u,c,lst,ns,n,TYPE)
  USE arnoldi_global,    ONLY    : iout
  IMPLICIT NONE
  INTEGER                                :: n, ns
  REAL*8, DIMENSION(n)                   :: phi
  REAL*8, DIMENSION(n,n)                 :: u
  REAL*8, DIMENSION(ns)                  :: c
  INTEGER,DIMENSION(ns)                  :: lst
  CHARACTER (LEN=*)                      :: TYPE
  INTEGER                                :: i
  phi=0.d0
  do i=1,ns
     phi(:) = phi(:) + u(:,lst(i))*c(i)
  END DO
  WRITE(iout,1) TYPE, ( lst(i),c(i), i=1,ns)
1    FORMAT(/,1X,'list of ',a1,' states and coefficients '//  &
    'superposed',(/,1X,5('state = ',i3,1X, 'coefficient = ',e15.5)))
END SUBROUTINE mk_phi

