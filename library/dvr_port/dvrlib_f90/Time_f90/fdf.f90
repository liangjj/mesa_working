!deck fdf.f
!***begin prologue     fdf
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            lobatto matrix elements of function
!***                   and first derivative.
!***
!***references

!***routines called
!***end prologue       fdf

  SUBROUTINE fdf(fa,dfb,wt,mat,n)
  IMPLICIT NONE
  INTEGER                            :: n
  INTEGER                            :: i, j
  REAL*8, DIMENSION(n,n)             :: fa, dfb, mat
  REAL*8, DIMENSION(n)               :: wt
  mat=0
  DO  i=1,n
    DO  j=1,n
      mat(i,j) = fa(i,i)*wt(i)*dfb(i,j)
    END DO
  END DO
  RETURN
END SUBROUTINE fdf



