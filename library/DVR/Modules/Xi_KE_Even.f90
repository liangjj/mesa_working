!deck Xi_KE_Even.f
!***begin prologue     Xi_KE_Even
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            angular KE elements for prolate coordinates.
!*** description       Here we consider the matrix elements
!***                   < f_i | d/dx (x*x - 1. ) d/dx | f_j>
!                      = - <  | [ df_i/dx df_j/dx ] ( 1 - x*x ) >
!                      where the last line comes from integrating by parts and dropping the
!                      vanishing surface term.       
!***references

!***routines called
!***end prologue       Xi_KE_Even

  SUBROUTINE Xi_KE_Even(tr,q,wt,f,df,n)
  USE dvr_global ,  ONLY  : iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                              :: n
  REAL*8, DIMENSION(:)                 :: q
  REAL*8, DIMENSION(:)                 :: wt
  REAL*8, DIMENSION(:,:)               :: f
  REAL*8, DIMENSION(:,:)               :: df
  REAL*8, DIMENSION(:,:)               :: tr
  REAL*8, DIMENSION(:)                 :: fac
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: k
  CALL rzero(tr,n*n)
  DO  i = 1, n     
      DO j = 1, i
         DO k = 1, n
            tr(i,j) = tr(i,j) - fac(k) * wt(k) * df(k,i) * df(k,j) 
         END DO
         tr(j,i) = tr(i,j)
      END DO
  END DO
END SUBROUTINE Xi_KE_Even
