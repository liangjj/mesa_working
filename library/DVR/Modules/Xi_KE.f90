!deck Xi_KE.f
!***begin prologue     Xi_KE
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
!***end prologue       Xi_KE

  SUBROUTINE Xi_KE(tr,q,wt,f,df,n,region)
  USE dvr_global ,  ONLY  : iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                              :: n
  REAL*8, DIMENSION(:)                 :: q
  REAL*8, DIMENSION(:)                 :: wt
  REAL*8, DIMENSION(:,:)               :: f
  REAL*8, DIMENSION(:,:)               :: df
  REAL*8, DIMENSION(:,:)               :: tr
  REAL*8, DIMENSION(:), ALLOCATABLE    :: fac
  REAL*8                               :: one = 1.d0
  REAL*8                               :: two = 2.d0
  INTEGER                              :: region
  CHARACTER (LEN=80)                   :: title
  CHARACTER (LEN=3)                    :: itoc
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: k
  CALL rzero(tr,n*n)
  ALLOCATE ( fac(n) )
  DO  i = 1, n
      fac(i) = ( q(i) * q(i) - one )
  END DO
  ij = 0
  DO  i = 1, n     
      DO j = 1, i
         DO k = 1, n
            tr(i,j) = tr(i,j) - fac(k) * wt(k) * df(k,i) * df(k,j) 
         END DO
      END DO
  END DO
  IF( m_max >= 0) THEN
      DO  i = 1, n
          DO j = 1, i
             DO k = 1, n
                tr(i,j) = tr(i,j) -  fac(k) * fac(k) * wt(k) * df(k,i) * df(k,j)
             END DO
             tr(i,j) = tr(i,j)  - q(i) * fac(i) * f(i,i) * df(i,j)    &
                                - q(j) * fac(j) * f(j,j) * df(j,i) )  
             tr(j,i) = tr(i,j)
          END DO
          tr(i,i) = tr(i,i)  - q(i) * q(i) * f(i,i) * f(i,i)
      END DO
  END IF
END SUBROUTINE Xi_KE
