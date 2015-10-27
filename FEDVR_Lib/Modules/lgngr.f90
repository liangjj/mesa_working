!*deck lgngr
      SUBROUTINE LGNGR(p,dp,ddp,x,y,nx,ny,prnt,drctv,type) 
!***begin prologue     lgngr
!***date written       940504   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G% 
!***purpose            lagrange polynomials at arbitrary points.
!***description
!***            
!               
!               
!***references
!
!***routines called
!
!***end prologue       lgngr
!
      IMPLICIT NONE
      REAL*8, DIMENSION (ny,nx)    :: p
      REAL*8, DIMENSION (ny,nx)    :: dp
      REAL*8, DIMENSION (ny,nx)    :: ddp
      REAL*8, DIMENSION(nx)        :: x
      REAL*8, DIMENSION(ny)        :: y
      REAL*8                       :: sn
      REAL*8                       :: ssn
      REAL*8                       :: fac
      LOGICAL                      :: prnt
      CHARACTER (LEN = 80)         :: title 
      CHARACTER (LEN = *)          :: drctv
      CHARACTER (LEN = *)          :: type
      INTEGER                      :: nx
      INTEGER                      :: ny
      INTEGER                      :: i
      INTEGER                      :: j
      INTEGER                      :: k
      INTEGER                      :: first
      INTEGER                      :: second
      INTEGER                      :: zerfac
      INTEGER                      :: inp
      INTEGER                      :: iout
      common /io/ inp, iout

!
!     generate polynomials and derivatives with respect to x
!
      p(:,:) = 1.d0
      IF (type /= 'even') THEN 
          DO i=1,ny
             zerfac = 0
             DO j=1,nx
                fac =  y(i) - x(j) 
                IF(abs(fac) <= 1.d-10) THEN
                   zerfac = j
                ENDIF  
             END DO
             DO j=1,nx
                DO k = 1, j-1
                   p(i,j) = p(i,j) * ( y(i) - x(k) )   &
                                   / ( x(j) - x(k) )
                END DO
                DO k=j+1,nx
                   p(i,j) = p(i,j) * ( y(i) - x(k) )   &
                                   / ( x(j) - x(k) )
                END DO
                IF( drctv /= 'functions-only') THEN
                    IF ( abs(p(i,j)) > 1.d-10) THEN
                         sn = 0.d0
                         ssn = 0.d0
                         DO k=1,j-1
                            fac = 1.d0/( y(i) - x(k) )
                            sn = sn + fac
                            ssn = ssn + fac*fac
                         END DO
                         DO k=j+1,nx
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
                         DO k=1,first-1
                            fac = 1.d0/( x(j) - x(k) )
                            sn = sn*fac*( y(i) - x(k) )
                            ssn = ssn + 1.d0/(y(i) - x(k))
                         END DO
                         DO k=first+1,second-1
                            fac = 1.d0/( x(j) - x(k) )
                            sn = sn*fac*( y(i) - x(k) )
                            ssn = ssn + 1.d0/( y(i) - x(k) )             
                         END DO
                         DO k=second+1,nx
                            fac = 1.d0/( x(j) - x(k) )
                            sn = sn*fac*( y(i) - x(k) )
                            ssn = ssn + 1.d0/( y(i) - x(k) )             
                         END DO
                         dp(i,j) = sn/( x(j) - x(zerfac) )
                         ddp(i,j) = 2.d0*ssn*dp(i,j)
                    END IF                    
                END IF
             END DO
!
          END DO
      ELSE 
          DO i=1,ny
             zerfac = 0
             DO j=1,nx
                fac =  y(i)*y(i) - x(j)*x(j) 
                IF(abs(fac) <= 1.d-10) THEN
                   zerfac = j
                ENDIF  
             END DO
             DO j=1,nx
                DO k=1,j-1
                   p(i,j) = p(i,j) * ( y(i)*y(i) - x(k)*x(k) )   &
                                / ( x(j)*x(j) - x(k)*x(k) )
                END DO
                DO k=j+1,nx
                   p(i,j) = p(i,j) * ( y(i)*y(i) - x(k)*x(k) )  &
                                / ( x(j)*x(j) - x(k)*x(k) )
                END DO
                IF( drctv /= 'functions-only') THEN
                    IF ( abs(p(i,j)) > 1.d-10) THEN
                         sn = 0.d0
                         ssn = 0.d0
                         DO k=1,j-1
                            fac = 1.d0/( y(i)*y(i) - x(k)*x(k) )
                            sn = sn + fac
                            ssn = ssn + fac*fac
                         END DO
                         DO k=j+1,nx
                            fac = 1.d0/( y(i)*y(i) - x(k)*x(k) )
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
                         DO k=1,first-1
                            fac = 1.d0/( x(j)*x(j) - x(k)*x(k) )
                            sn = sn*fac*( y(i)*y(i) - x(k)*x(k) )
                            ssn = ssn + 1.d0/(y(i)*y(i) - x(k)*x(k))
                         END DO
                         DO k=first+1,second-1
                            fac = 1.d0/( x(j)*x(j) - x(k)*x(k) )
                            sn = sn*fac*( y(i)*y(i) - x(k)*x(k) )
                            ssn = ssn + 1.d0/( y(i)*y(i) - x(k)*x(k) )             
                         END DO
                         DO k=second+1,nx
                            fac = 1.d0/( x(j)*x(j) - x(k)*x(k) )
                            sn = sn*fac*( y(i)*y(i) - x(k)*x(k) )
                            ssn = ssn + 1.d0/( y(i)*y(i) - x(k)*x(k) )             
                         END DO
                         dp(i,j) = sn/( x(j)*x(j) - x(zerfac)*x(zerfac) )
                         ddp(i,j) = 2.d0*ssn*dp(i,j)
                         ddp(i,j) = 2.d0*dp(i,j) + 4.d0 * y(i) * y(i) * ddp(i,j) 
                         dp(i,j) = 2.d0 * y(i) * dp(i,j)
                    END IF                    
                END IF
             END DO
          END DO
!
      END IF
      IF(prnt) THEN
          title='polynomials'
          call prntfm(title,p,ny,nx,ny,nx,iout)
          IF(drctv.ne.'functions-only') then
             title='derivative of polynomials'
             call prntfm(title,dp,ny,nx,ny,nx,iout)
             title='second derivative of polynomials'
             call prntfm(title,ddp,ny,nx,ny,nx,iout)
          END IF
      END IF
      END SUBROUTINE Lgngr
