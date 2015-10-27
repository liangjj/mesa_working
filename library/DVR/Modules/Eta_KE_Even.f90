!deck KE_Mat_Even
!**begin prologue      KE_Mat_Even
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!*** description              
!***references

!***routines called
!***end prologue       KE_Mat_Even

  SUBROUTINE KE_Mat_Even(tr,df,q,wt,fac,n)
  IMPLICIT NONE
  INTEGER                              :: n
  REAL*8, DIMENSION(:,:)               :: tr
  REAL*8, DIMENSION(:,:)               :: df
  REAL*8, DIMENSION(:)                 :: q
  REAL*8, DIMENSION(:)                 :: wt
  REAL*8, DIMENSION(:)                 :: fac


  REAL*8, DIMENSION(:)                 :: fac
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: k
  CALL rzero(tr,n*n)
  END DO
  DO i = 1, n
     DO j = 1, i
        DO k = 1, n
           tr(i,j) = tr(i,j) - fac(k) * wt(k) * df(k,i) * df(k,j) 
        END DO
        tr(j,i) = tr(i,j)
     END DO
  END DO
END SUBROUTINE Eta_KE_Even
