!deck vexp123.f
!***begin prologue     vexp123
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            additive pair potential
!***
!***references
!***routines called
!***end prologue       vexp123
  SUBROUTINE vexp123(v,q1,q2,q3,a12,a13,a23, b12,b13,b23,n1,n2,n3,prn)
  USE input_output
  INTEGER                                :: n1, n2, n3
  REAL*8, DIMENSION(n3,n2,n1)            :: v
  REAL*8, DIMENSION(n1)                  :: q1
  REAL*8, DIMENSION(n2)                  :: q2
  REAL*8, DIMENSION(n3)                  :: q3
  REAL*8                                 :: a12, a13, a23, b12, b13, b23
  LOGICAL                                :: prn
  REAL*8                                 :: add
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=16)                     :: fptoc
  DO  i=1,n1
      DO  j=1,n2
          add = a12 * EXP ( - b12 * ABS( q2(j) - q1(i) ) )
          v(:,j,i) = v(:,j,i) + add
      END DO
  END DO
  DO  i=1,n1
      DO  k=1,n3
          add = a13 * EXP ( - b13 * ABS( q3(k) - q1(i) ) )
          v(k,:,i) = v(k,:,i) + add
      END DO
  END DO
  DO  j=1,n2
      DO  k=1,n3
         add = a23 * EXP ( - b23 * ABS( q3(k) - q2(j) ) )
         v(k,j,:) = v(k,j,:) + add
      END DO
  END DO
  IF(prn) THEN
     DO  i=1,n1
         title='sum of two dimensional exponentials x= ' //fptoc(q1(i))
         CALL prntrm(title,v(1,1,i),n3,n2,n3,n2,iout)
     END DO
  END IF
END SUBROUTINE vexp123



