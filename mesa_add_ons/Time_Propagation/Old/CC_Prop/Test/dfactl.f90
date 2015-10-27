!deck dfactl
!begin prologue     dfactl
!date written       920405   (yymmdd)
!revision date      yymmdd   (yymmdd)
!keywords           double factorials
!author             schneider, barry (nsf)
!source             mylib
!purpose            calculate double factorials from zero to n
!description        
!                   
!                   
!references         
!routines called    
!end prologue      factl
      subroutine dfactl(dfact,n)
      IMPLICIT NONE
      INTEGER                          :: n
      INTEGER                          :: i
      REAL*8, DIMENSION(0:n)           :: dfact
      dfact(0)=1.d0
      dfact(1)=1.d0
      dfact(2)=3.d0
      DO i=3,n
         dfact(i) = ( i + i - 1) * dfact(i-1) 
      END DO
  END SUBROUTINE dfactl
