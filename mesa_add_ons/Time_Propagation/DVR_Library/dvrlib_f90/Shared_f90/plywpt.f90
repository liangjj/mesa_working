!deck plywpt.f
!***begin prologue     plywpt
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            lobatto points, weights, polynomials
!***                   and their first and second derivatives.
!***
!***references
!***routines called
!***end prologue       plywpt
  SUBROUTINE plywpt(q,wt,p,dp,ddp,cordsys,qdtyp,edge,n,nq,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n, nq 
  REAL*8, DIMENSION(n)                   :: q, wt
  REAL*8, DIMENSION(n,n)                 :: p, dp, ddp
  CHARACTER (LEN=*)                      :: cordsys, qdtyp
  REAL*8, DIMENSION(2)                   :: edge
  LOGICAL                                :: prn(*)
  CHARACTER (LEN=80)                     :: title
  CALL drvply(q,wt,p,dp,ddp,edge,cordsys,qdtyp,n-1,nq,prn)
  IF(prn(1)) THEN
     title='points'
     CALL prntrm(title,q,n,1,n,1,iout)
     title='weights'
     CALL prntrm(title,wt,n,1,n,1,iout)
  END IF
  IF(prn(2)) THEN
     title='polynomials'
     CALL prntrm(title,p,n,n,n,n,iout)
     title='first derivative polynomials'
     CALL prntrm(title,dp,n,n,n,n,iout)
     title='second derivative polynomials'
     CALL prntrm(title,ddp,n,n,n,n,iout)
  END IF
END SUBROUTINE plywpt


