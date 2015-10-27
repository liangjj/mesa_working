!deck tplmat.f
!
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            lobatto matrix elements.
!***
!***references

!***routines called
!***end prologue       tplmat

  SUBROUTINE tplmat(p,dp,wt,mat,n)
  IMPLICIT NONE
  INTEGER n
  REAL*8, DIMENSION(n,n)                :: p, dp, mat
  REAL*8, DIMENSION(n)                  :: wt
  CHARACTER (LEN=80) :: title
  CALL fdf(p,dp,wt,mat,n)
  RETURN
END SUBROUTINE tplmat






