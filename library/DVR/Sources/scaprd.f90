!deck scaprd
!***begin prologue     scaprd
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           
!***author             schneider, barry (nsf)
!***source
!***purpose
!***
!***
!***description
!***
!***
!***
!***references

!***routines called
!***end prologue       scaprds
  Function Scaprd (va,vb,wt,rwt,n)
  IMPLICIT NONE
  INTEGER                 :: n
  INTEGER                 :: i
  REAL*8, DIMENSION(n)    :: va
  REAL*8, DIMENSION(n)    :: vb
  REAL*8                  :: wt
  REAL*8                  :: rwt
  REAL*8                  :: scaprd
  scaprd=0.d0
  DO  i=1,n
      scaprd=scaprd+va(i)*wt(i)*rwt(i)*vb(i)
  END DO
END FUNCTION Scaprd
