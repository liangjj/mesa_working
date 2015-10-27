!deck legend
!begin prologue     legend
!date written       880721   (yymmdd)
!revision date      yymmdd   (yymmdd)
!keywords           legend, link 2702, legendre functions
!author             schneider, barry (lanl)
!source             m2702
!purpose            legendre functions
!description        calculation of p(l,m) functions
!references         none
!                      plm are the legendre functions l=m to l=lmax    
!                      x are the values of cos(theta)
!                      dfct are the factorials from 0 to maxfac
!                      ddfct are the double factorials from 0 to maxfac    
!routines called
!end prologue       legend
      subroutine legend (plm,x,dfct,ddfct,npt,lmax,m,maxfac)
      IMPLICIT NONE
      INTEGER                          :: npt, lmax, m, maxfac
      REAL*8, DIMENSION(npt)           :: x
      REAL*8, DIMENSION(0:maxfac)      :: dfct, ddfct
      REAL*8, DIMENSION(npt,m:maxfac)  :: plm
      REAL*8                           :: fm, facx, f1
      INTEGER                          :: i, n1, n2, n3
!----------------------------------------------------------------------
!           start recursion with plm(m,m) and plm(m+1,m)               
!                      and recur upward                                
!----------------------------------------------------------------------
      plm(1:npt,m:lmax) = 0.d0
!
!     Get first plm
!
      IF ( m == 0) THEN
          plm(:,m) = ddfct(m)
      ELSE
          fm=.5d+00*m
          plm(:,m) = ddfct(m) * (1.d+00-x(:)*x(:))**fm
      END IF
!
!     Get second plm to do the recursion
!
      IF (lmax /= m) THEN
          plm(:,m+1) = ( m + m + 1) * x(:) * plm(:,m)
!
!         Recur if required
!
          IF (lmax /= m+1) THEN
              n1=2
              n2=m+m+3
              n3=n2-2
              DO i=m+2,lmax
                 plm(:,i) = ( n2 * x(:) * plm(:,i-1) - n3 * plm(:,i-2) ) / n1
                 n1=n1+1
                 n2=n2+2
                 n3=n3+1
             END DO
          END IF
      END IF
  END SUBROUTINE legend
