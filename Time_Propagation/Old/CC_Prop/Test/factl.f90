!*deck factl
!***begin prologue     factl
!***date written       920405   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           factorials
!***author             schneider, barry (nsf)
!***source             mylib
!***purpose            calculate factorials from zero to n
!***description        
!***                   
!***                   
!
!***references         
!***routines called    
!***end prologue      factl
      subroutine factl(fact,n)
      USE input_output
      IMPLICIT NONE
      INTEGER                             :: n, i
      real*8, DIMENSION(0:n)              :: fact
      INTEGER                             :: maxn = 100
      IF(n > maxn) THEN
         call lnkerr('factorial will overflow')
      END IF
      fact(0)=1.d0
      DO i=1,n
         fact(i) = i * fact(i-1) 
      END DO
  END SUBROUTINE factl
