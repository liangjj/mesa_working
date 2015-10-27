!deck todelta.f
!***begin prologue     todelta
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            renormalization
!***
!***description
!***references
!***routines called
!***end prologue       todelta
  SUBROUTINE todelta(p,dp,ddp,wt,n)
  IMPLICIT NONE
  INTEGER                                :: n, i  
  REAL*8, DIMENSION(n,n)                 :: p, dp, ddp
  REAL*8, DIMENSION(n)                   :: wt
  REAL*8                                 :: norm
  DO  i=1,n
      norm=SQRT(wt(i))
      p(:,i)= p(:,i)*norm
      dp(:,i)= dp(:,i)*norm
      ddp(:,i)= ddp(:,i)*norm
  END DO
END SUBROUTINE todelta
