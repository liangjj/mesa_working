!deck vwel123.f
!***begin prologue     vwel123
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            two dimensional potential well
!***
!***references
!***routines called
!***end prologue       vwel123
  SUBROUTINE vwel123(v,q1,q2,q3,a12,a13,a23, d12,d13,d23,n1,n2,n3,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3
  REAL*8, DIMENSION(n3,n2,n1)            :: v
  REAL*8, DIMENSION(n1)                  :: q1
  REAL*8, DIMENSION(n2)                  :: q2
  REAL*8, DIMENSION(n3)                  :: q3
  REAL*8                                 :: a12, a13, a23, d12, d13, d23
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=16)                     :: fptoc
  INTEGER                                :: i, j, k
  DO  i=1,n1
      DO  j=1,n2
          IF( ABS(q1(i)-q2(j)) <= a12) THEN
             v(:,j,i) = v(:,j,i) + d12
          END IF
      END DO
  END DO
  DO  i=1,n1
      DO  k=1,n3
          IF( ABS(q1(i)-q3(k)) <= a13) THEN
              v(k,:,i) = v(k,:,i) + d13
          END IF
      END DO
  END DO
  DO  j=1,n2
      DO  k=1,n3
          IF( ABS(q2(j)-q3(k)) <= a23) THEN
             v(k,j,:) = v(k,j,:) + d23
          END IF
      END DO
  END DO
  IF(prn) THEN
     DO  i=1,n1
         title='sum of two dimensional wells x= '//fptoc(q1(i))
         CALL prntrm(title,v(1,1,i),n3,n2,n3,n2,output)
     END DO
  END IF
1    FORMAT(/,9X,'coordinate',8X,'potential')
2    FORMAT(5X,e15.8,3X,e15.8)
END SUBROUTINE vwel123



