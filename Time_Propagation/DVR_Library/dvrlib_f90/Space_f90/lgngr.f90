!deck lgngr
  SUBROUTINE lgngr(p,dp,ddp,x,y,nx,ny,prnt,drctv)
!***begin prologue     lgrply
!***date written       940504   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            lagrange polynomials at arbitrary points.
!***description
!***
!***references
!***routines called
!***end prologue       lgngr
  USE input_output 
  IMPLICIT NONE
  INTEGER                                :: nx, ny
  REAL*8, DIMENSION(ny,nx)               :: p, dp, ddp
  REAL*8, DIMENSION(nx)                  :: x
  REAL*8, DIMENSION(ny)                  :: y
  LOGICAL                                :: prnt
  CHARACTER (LEN=*)                      :: drctv
  REAL*8                                 :: sn, ssn, fac
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: i, j, k, zerfac, first, second
!
!     generate polynomials and derivatives with respect to x
!
  DO  i=1,ny
      zerfac = 0
      DO  j=1,nx
          fac =  y(i) - x(j)
          IF(ABS(fac) <= 1.d-10) THEN
             zerfac = j
          END IF
      END DO
      DO  j=1,nx
          p(i,j) = 1.d0
          DO  k=1,j-1
             p(i,j) = p(i,j)*( y(i) - x(k) ) / ( x(j) - x(k) )
          END DO
          DO  k=j+1,nx
              p(i,j) = p(i,j)*( y(i) - x(k) ) / ( x(j) - x(k) )
          END DO
          IF(drctv /= 'functions-only') THEN
             IF(ABS(p(i,j)) > 1.d-10) THEN
                sn = 0.d0
                ssn = 0.d0
                DO  k=1,j-1
                    fac = 1.d0/( y(i) - x(k) )
                    sn = sn + fac
                    ssn = ssn + fac*fac
               END DO
               DO  k=j+1,nx
                   fac = 1.d0/( y(i) - x(k) )
                   sn = sn + fac
                   ssn = ssn + fac*fac
               END DO
               dp(i,j) = sn*p(i,j)
               ddp(i,j) = sn*dp(i,j) - ssn*p(i,j)
             ELSE
               first=j
               second=zerfac
               IF(first > second) THEN
                  first=zerfac
                  second=j
               END IF
               sn = 1.d0
               ssn = 0.d0
               DO  k=1,first-1
                  fac = 1.d0/( x(j) - x(k) )
                  sn = sn*fac*( y(i) - x(k) )
                  ssn = ssn + 1.d0/(y(i) - x(k))
               END DO
               DO  k=first+1,second-1
                   fac = 1.d0/( x(j) - x(k) )
                   sn = sn*fac*( y(i) - x(k) )
                   ssn = ssn + 1.d0/( y(i) - x(k) )
               END DO
               DO  k=second+1,nx
                   fac = 1.d0/( x(j) - x(k) )
                   sn = sn*fac*( y(i) - x(k) )
                   ssn = ssn + 1.d0/( y(i) - x(k) )
               END DO
               dp(i,j) = sn/( x(j) - x(zerfac) )
               ddp(i,j) = 2.d0*ssn*dp(i,j)
             END IF
          END IF
      END DO
  END DO
  IF(prnt) THEN
     title='polynomials'
     CALL prntfm(title,p,ny,nx,ny,nx,iout)
     IF(drctv /= 'functions-only') THEN
        title='derivative of polynomials'
        CALL prntfm(title,dp,ny,nx,ny,nx,iout)
        title='second derivative of polynomials'
        CALL prntfm(title,ddp,ny,nx,ny,nx,iout)
     END IF
  END IF
END SUBROUTINE lgngr















