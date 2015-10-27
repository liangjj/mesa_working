!deck fourier.f
!***begin prologue     fourier
!***date written       000702   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr
!***
!***author             schneider, b. i.(nsf)
!***source             fourier
!***purpose            Points, and weights for fourier dvr
!***description        
!***                   
!***                   
!***                   
!***                   

!***references         

!***routines called    iosys, util and mdutil
!***end prologue       fourier
  SUBROUTINE fourier(q,wt,f,edge,nord)
  USE dvr_global,       ONLY  : inp, iout, box, deltax, pi
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nord, i, j, k
  REAL*8, DIMENSION(nord)                :: q, wt
  REAL*8, DIMENSION(nord,nord)           :: f
  REAL*8, DIMENSION(2)                   :: edge
  REAL*8                                 :: fac
  CHARACTER(LEN=3)                       :: itoc
!
  box=edge(2)-edge(1)
  deltax=box/nord
  wt = deltax
  j=(nord-1)/2
  fac = edge(1) + .5d0*box
  k=-j
  DO i=1,nord
     q(i) = deltax * k + fac   
     k=k+1
  END DO
  f=0.d0
  fac=1.d0/sqrt(deltax)
  DO i=1,nord
     f(i,i) = fac
  END DO
!
END SUBROUTINE fourier
