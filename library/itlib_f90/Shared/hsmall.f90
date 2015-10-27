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

  SUBROUTINE hsmall(begin,END,action,prnt)
  USE io
  USE dvd_prnt
  USE dvd_global
  IMPLICIT NONE
  INTEGER                                :: begin
  INTEGER                                :: END
  CHARACTER (LEN=*)                      :: action
  LOGICAL                                :: prnt
  REAL*8                                 :: sdot
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: i, j
  IF(action == 'initialize') THEN
     DO  i=1,END
         DO  j=1,i
             b(i,j) = sdot(n,vec(1,i),1,hvec(1,j),1)
             b(j,i) = b(i,j)
             btmp(i,j) = b(i,j)
             btmp(j,i) = b(j,i)
         END DO
     END DO
  ELSE IF(action == 'fill') THEN
     DO  i=1,END
         DO  j=begin,END
             b(i,j) = sdot(n,vec(1,i),1,hvec(1,j),1)
             b(j,i) = b(i,j)
         END DO
     END DO
     DO  i=1,END
         DO  j=1,i
             btmp(i,j) = b(i,j)
             btmp(j,i) = b(i,j)
         END DO
     END DO
   END IF
   IF(prnt) THEN
      title='small matrix'
      CALL prntfm(title,btmp,END,maxvec,END,maxvec,iout)
   END IF
END SUBROUTINE hsmall

