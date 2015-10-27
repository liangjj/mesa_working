!***********************************************************************
! Gauss_Quadrature
!**begin prologue     Gauss_Quadrature
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           FEDVR, grid, polynomials
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Calculate the regional FEDVR points, weights,
!***                  polynomials and normalizations factors as well as
!***                  some global poibnts and weights.
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***end prologue      Gauss_Quadrature
!***********************************************************************
!***********************************************************************
                           MODULE Gauss_Quadrature
                           USE Grid_Defined_Types
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck gauss.f
!***begin prologue     gauss
!***date written       000702   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr
!***
!***author             schneider, b. i.(nsf)
!***source             gauss
!***purpose            Points, weights and coordinate functions for generalized
!***                   Gauss quadratures.
!***description        Lanczos recursion using a reference weight function
!***                   is used to generate the points and weights of Gauss quadratures
!***                   for generalized weight functions.  The eigenvector matrix
!***                   of the tridiagonal matrix is used to compute the
!***                   coordinate functions, their first and second derivatives.

!***references         see papers and notes appended.

!***routines called    iosys, util and mdutil
!***end prologue       gauss

!     This is the main library routine to compute the orthogonal and
!     coordinate DVR functions for general weight functions.  The approach
!     is to use a reference quadrature to compute the $\alpha$ and $\beta$
!     recursion coefficients satisfied by the orthogonal polynonials.  The
!     three term recursion relationship is then diagonalized to obtain the
!     generalized points and weights.  The eigenvectors are used to transform
!     from the orthogonal polynomial to coordinate representation.  First
!     and second derivatives are also computed by taking the derivative of the
!     recursion relationship.

  SUBROUTINE gauss(q,wt,edge,p,dp,ddp,type_quadrature,fixed_point,n,print)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:)                :: q
  REAL(idp), DIMENSION(:)                :: wt
  REAL(idp), DIMENSION(:)                :: edge
  REAL(idp), OPTIONAL, DIMENSION(:,:)    :: p
  REAL(idp), OPTIONAL, DIMENSION(:,:)    :: dp
  REAL(idp), OPTIONAL, DIMENSION(:,:)    :: ddp
  REAL(idp),  DIMENSION(:), ALLOCATABLE  :: b
  CHARACTER(LEN=*)                       :: type_quadrature
  INTEGER                                :: fixed_point
  INTEGER                                :: n
  REAL(idp), DIMENSION(2)                :: endpts
  REAL(idp), DIMENSION(2)                :: ptfix
  LOGICAL, OPTIONAL                      :: print
  DATA ptfix / -1.d0, 1.d0 /
!
!
! If the weight function is a one of the classical weight functions the
! points and weights are known analytically and after computing them we
! go directly to getting the coordinate functions.
!
  ALLOCATE( b(1:n) )
  endpts(:)=edge(:)
  IF (type_quadrature == "gauss") THEN
      CALL gaussq('one',n,0.d0,0.d0,0,ptfix,b,q,wt)
  ELSE IF ( type_quadrature == "radau") THEN
      IF (fixed_point == 1) THEN
          ptfix(1) = -1.d0
      ELSE IF ( fixed_point == 2) THEN
         ptfix(1) = 1.d0
      END IF
      CALL gaussq('one',n,0.d0,0.d0,1,ptfix,b,q,wt)
  ELSE IF ( type_quadrature == "lobatto" ) THEN
      ptfix(1) = -1.d0
      ptfix(2) = 1.d0
      CALL gaussq('one',n,0.d0,0.d0,2,ptfix,b,q,wt)
  END IF
  CALL cnvtpt(q,wt,edge,n)
  IF(PRESENT(print) == .true.) THEN
     call Print_Matrix(type_real_vector,q,title='Final Nodes from Gauss')
     call Print_Matrix(type_real_vector,wt,title='Final Weights from Gauss')
  END IF
  DEALLOCATE(b)
!  
!  
! The DVR library assumes that the polynomials are $\delta$
! functions at the quadrature points.  Convert to this normalization
!
!
END SUBROUTINE gauss
!***********************************************************************
!***********************************************************************
!*deck lgngr
   SUBROUTINE LGNGR(p,dp,ddp,x,y,nx,ny,type,drctv,print) 
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
  REAL*8, DIMENSION (ny,nx)           :: p
  REAL*8, DIMENSION (ny,nx)           :: dp
  REAL*8, DIMENSION (ny,nx)           :: ddp
  REAL*8, DIMENSION(nx)               :: x
  REAL*8, DIMENSION(ny)               :: y
  REAL*8, DIMENSION(:), ALLOCATABLE   :: xt
  REAL*8, DIMENSION(:), ALLOCATABLE   :: yt
  REAL*8                              :: sn
  REAL*8                              :: ssn
  REAL*8                              :: fac
  LOGICAL, OPTIONAL                   :: print
  CHARACTER (LEN = 80)                :: title 
  CHARACTER (LEN = *), OPTIONAL       :: drctv
  CHARACTER (LEN = *), OPTIONAL       :: type
  INTEGER                             :: nx
  INTEGER                             :: ny
  INTEGER                             :: i
  INTEGER                             :: j
  INTEGER                             :: k
  INTEGER                             :: first
  INTEGER                             :: second
  INTEGER                             :: zerfac
  INTEGER                             :: inp
  INTEGER                             :: iout
!
!     generate polynomials and derivatives with respect to x
!
  p(:,:) = 1.d0
  IF (present(type) ) THEN 
      ALLOCATE(xt(nx),yt(ny))
      xt(:) = x(:)
      yt(:) = y(:)
      x(:) = x(:) * x(:)
      y(:) = y(:) * y(:)
  END IF
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
        IF(present(drctv) ) THEN
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
!
     END DO
  END DO
!
  IF (present(type)) THEN 
      DO i=1,ny
         ddp(i,:) = 2.d0*dp(i,:) + 4.d0 * yt(i) * yt(i) * ddp(i,:) 
         dp(i,:) = 2.d0 * yt(i) * dp(i,:)
      END DO
      x(:) = xt(:)
      y(:) = yt(:)
      DEALLOCATE(xt,yt)
!
  END IF
  IF(present(print)) THEN
     call Print_Matrix(type_real_matrix,p,ny,nx,title='Polynomials')
     IF(present(drctv)) THEN
        call Print_Matrix(type_real_matrix,dp,ny,nx,title='First Derivative of Polynomials')
        call Print_Matrix(type_real_matrix,ddp,ny,nx,title='Second Derivative of Polynomials')
     END IF
  END IF
  END SUBROUTINE Lgngr
!***********************************************************************
!***********************************************************************
!deck cnvtpt.f
  SUBROUTINE cnvtpt(pt,wt,endpts,n)
  IMPLICIT NONE
  INTEGER                                :: n
  REAL(idp), DIMENSION(n)                :: pt
  REAL(idp), DIMENSION(n)                :: wt
  REAL(idp), DIMENSION(2)                :: endpts
  REAL(idp)                              :: f1
  REAL(idp)                              :: f2
  f1 = ( endpts(2)-endpts(1) )*.5D0
  f2 = ( endpts(1) + endpts(2) )*.5D0
  pt =  f1*pt + f2
  wt = wt*f1
END SUBROUTINE cnvtpt
!**********************************************************************
!**********************************************************************
!deck cpoly.f
!***begin prologue     cpoly
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate coordinate functions.
!***
!***description
!***references
!***routines called
!***end prologue       cpoly

  SUBROUTINE cpoly(cp,dcp,ddcp,pt,n,prn)
  IMPLICIT NONE
  INTEGER                             :: n
  REAL(idp), DIMENSION(:,:)           :: cp
  REAL(idp), DIMENSION(:,:)           :: dcp
  REAL(idp), DIMENSION(:,:)           :: ddcp
  REAL(idp), DIMENSION(:)             :: pt
  LOGICAL, OPTIONAL                   :: prn
  CHARACTER (LEN=80)                  :: title
  CALL lgngr(cp,dcp,ddcp,pt,pt,n,n,drctv='on')
  IF(PRESENT(prn) == .true.) THEN
     call Print_Matrix(type_real_matrix,cp,n,n,title='Coordinate Function')
     call Print_Matrix(type_real_matrix,dcp,n,n,title='First Derivative of Coordinate Function')
     call Print_Matrix(type_real_matrix,ddcp,n,n,title='Second Derivative of Coordinate Function')
   END IF
END SUBROUTINE cpoly
!***********************************************************************
!***********************************************************************
          END MODULE Gauss_Quadrature
!***********************************************************************
!***********************************************************************
