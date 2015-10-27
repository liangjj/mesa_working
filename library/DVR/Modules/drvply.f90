!deck drvply.f
!***begin prologue     drvply
!***date written       000702   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, gauss-qudrature, polynomials
!***
!***author             schneider, b. i.(nsf)
!***source             drvply
!***purpose            Simplified version of library routine drvply.
!***                   
!***description        Lanczos recursion using a reference weight function
!***                   is used to generate the points and weights of Gauss quadratures
!***                   for generalized weight functions.  The eigenvector matrix
!***                   of the tridiagonal matrix is used to compute the
!***                   coordinate functions, their first and second derivatives.

!***references         see papers and notes appended.

!***routines called    iosys, util and mdutil
!***end prologue       drvply

!     This is the main library routine to compute the orthogonal and
!     coordinate DVR functions.
!
 SUBROUTINE drvply(q,wt,p,dp,ddp,edge,typwt,nord,reg_number)
  USE dvr_global,       ONLY  : inp, iout, l_val, m_val, nfix, nreg, fix
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nord
  REAL*8, DIMENSION(nord)                :: q, wt
  REAL*8, DIMENSION(nord,nord)           :: p, dp, ddp
  REAL*8, DIMENSION(2)                   :: edge
  CHARACTER (LEN=*)                      :: typwt
  CHARACTER (LEN=80)                     :: title
  REAL*8, DIMENSION(2)                   :: endpts, ptse
  REAL*8, DIMENSION(2)                   :: ptfix
  REAL*8, DIMENSION(:), ALLOCATABLE      :: a, b
  DATA ptfix / -1.d0, 1.d0 /
!
!
!
  ALLOCATE( a(nord), b(nord) )
  endpts=edge
  ptse = endpts
  IF (nreg == 1) THEN
      n_fixed=0
      IF(nfix == 1) THEN
         n_fixed=1       
         ptfix(1) = -1.d0
         ptse(1) = endpts(1)
         IF(fix(2)) THEN
            ptfix(1) = 1.d0
            ptse(1) = endpts(2)
         END IF
      ELSE IF(nfix == 2) THEN
         n_fixed=2      
         ptfix(1) = -1.d0
         ptfix(2) = 1.d0
         ptse = endpts
      END IF
  ELSE
      IF(reg_number == 1) THEN 
         n_fixed=1       
         ptfix(1)=1.d0
         ptse(1) = endpts(2)
         IF(fix(1)) THEN
            n_fixed=2
            ptfix(1)=-1.d0
            ptfix(2)=1.d0
            ptse = endpts
         END IF
      ELSE IF(reg_number == nreg) THEN 
            n_fixed=1   
            ptfix(1)=-1.d0  
            ptse(1) = endpts(1)  
            IF(fix(2)) THEN
               n_fixed=2
               ptfix(1)=-1.d0
               ptfix(2)=1.d0
               ptse = endpts
            END IF
      ELSE
         n_fixed=2
         ptfix(1)=-1.d0
         ptfix(2)=1.d0
         ptse = endpts
      END IF
  END IF
!
! For finite regions
!
  write(iout,1) typwt
  CALL gaussq(typwt,nord,0.d0,0.d0,n_fixed,ptfix,b,q,wt)
  CALL cnvtpt(q,wt,edge,nord)
!
  CALL cpoly(p,dp,ddp,q,a,nord-1,nord,prn(5))
!  
! The DVR library assumes that the polynomials are $\delta$
! functions at the quadrature points.  Convert to this normalization
!
  DEALLOCATE(a,b)  
!
1    FORMAT(/,1X,'Generating Point and Weights for a Weight Function = ',a16)
END SUBROUTINE drvply
