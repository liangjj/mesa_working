!deck hsmall.f
!***begin prologue     hsmall
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           small davidson matrix
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       hsmall

  SUBROUTINE hsmall(prnt)
  USE io
  USE dvd_prnt
  USE dvd_global
  IMPLICIT NONE
  LOGICAL                                :: prnt
  REAL*8                                 :: sdot
  INTEGER                                :: i, j
  IF(drctv == 'initialize') THEN
     DO  i=1,size
         DO  j=1,i
             b(i,j) = sdot(n_dvd,vec(1,i),1,hvec(1,j),1)
             b(j,i) = b(i,j)
         END DO         
     END DO
     bwrk = b
  ELSE IF(drctv == 'fill') THEN
     DO  i=1,size
         DO  j=begin,size
             b(i,j) = sdot(n_dvd,vec(1,i),1,hvec(1,j),1)
             b(j,i) = b(i,j)
         END DO
     END DO
     DO  i=1,size
         DO  j=1,i
             bwrk(i,j) = b(i,j)
             bwrk(j,i) = b(i,j)
         END DO
     END DO
   END IF
   IF(prnt) THEN
      title='small matrix'
      CALL prntfm(title,bwrk,size,maxvec,size,maxvec,iout)
   END IF
END SUBROUTINE hsmall

