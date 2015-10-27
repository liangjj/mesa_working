!***********************************************************************
! CC_Prop_Module
!**begin prologue     CC_Prop_Module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Arnoldi, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the major subroutines to compute
!***                  the two electron integrals in a DVR product basis.  
!***description       
!***                  Two electron integrals in a product DVR basis are
!***                  computed by first solving the Poisson equation
!***                  for the density and then performing the second
!***                  integration by DVR quadrature.
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!**references
!**modules needed     See USE statements below
!**comments           In this portable version I have disabled all unnecessary
!**                   writing to files.  The original Fortran is commented out.
!**                   
!**end prologue       CC_Prop_Module
!***********************************************************************
!***********************************************************************
                           MODULE CC_Prop_Module
                           USE prop_prnt
                           USE dvrprop_global
                           USE dvr_shared
                           USE dvr_global
                        IMPLICIT NONE
!***********************************************************************
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!*deck CC_Prop_Driver
!***begin prologue     CC_Prop_Driver
!***date written       920405   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           two electron integrals
!***author             schneider, barry (nsf)
!***source             
!***purpose            driver for two electron integral calculation
!***description        
!***                   
!***                   
!
!***references         
!***routines called    
!***end prologue    CC_Prop_Driver
      subroutine CC_Prop_Driver
      IMPLICIT NONE
      REAL*8, DIMENSION(:),   ALLOCATABLE         :: fact, dfact
      REAL*8, DIMENSION(:,:), ALLOCATABLE         :: sin_m, cos_m
      REAL*8, DIMENSION(:,:), ALLOCATABLE         :: v_two
      TYPE legendre_functions
        REAL*8, DIMENSION(:,:),                             &
                POINTER                           :: plm
      END TYPE legendre_functions
      TYPE(legendre_functions),    DIMENSION(:),            &
                                   ALLOCATABLE    :: plm_mat
      INTEGER                                     :: i, n_tri
      ALLOCATE(fact(0:m_max), dfact(0:m_max),               &
               sin_m(nphy(3),m_max), cos_m(nphy(3),m_max))
      ALLOCATE(plm_mat(0:m_max))
      DO i=0,m_max
         ALLOCATE(plm_mat(i)%plm(nphy(2),i:l_max))
      END DO
!
!     Compute the needed P(l,m,cos(theta)) and Sin(m*phi)/Cos(m*phi) functions
!
      call factl(fact,m_max)
      call dfactl(dfact,m_max)
      DO i=0,m_max
         Call sin_cos(sin_m(:,i),cos_m(:,i),grid(3)%pt,nphy(3),i)
      END DO
      DO i=0,m_max
         Call legend(plm_mat(i)%plm,grid(2)%pt,dfact,nphy(2),l_max,i)
      END DO
      IF(.not.onel) THEN
         n_tri = nphy(1) * ( nphy(1) + 1 )/2
         ALLOCATE ( v_two(n_tri,0:l_max))      
         Call Two_Mat(v_two,grid(1)%ke,grid(1)%v,grid(1)%pt,grid(1)%wt,nphy(1))
      END IF
  END SUBROUTINE CC_Prop_Driver
!***********************************************************************
!*deck factl
!***begin prologue     factl
!***date written       920405   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           factorials
!***author             schneider, barry (nsf)
!***source             mylib
!***purpose            calculate factorials from zero to n
!***description        
!***                   
!***                   
!
!***references         
!***routines called    
!***end prologue      factl
      subroutine factl(fact,n)
      USE input_output
      IMPLICIT NONE
      INTEGER                             :: n, i
      real*8, DIMENSION(0:n)              :: fact
      INTEGER                             :: maxn = 100
      IF(n > maxn) THEN
         call lnkerr('factorial will overflow')
      END IF
      fact(0)=1.d0
      DO i=1,n
         fact(i) = i * fact(i-1) 
      END DO
  END SUBROUTINE factl
!***********************************************************************
!deck dfactl
!begin prologue     dfactl
!date written       920405   (yymmdd)
!revision date      yymmdd   (yymmdd)
!keywords           double factorials
!author             schneider, barry (nsf)
!source             mylib
!purpose            calculate double factorials from zero to n
!description        
!                   
!                   
!references         
!routines called    
!end prologue      factl
      subroutine dfactl(dfact,n)
      IMPLICIT NONE
      INTEGER                          :: n
      INTEGER                          :: i
      REAL*8, DIMENSION(0:n)           :: dfact
      dfact(0)=1.d0
      dfact(1)=1.d0
      dfact(2)=3.d0
      DO i=3,n
         dfact(i) = ( i + i - 1) * dfact(i-1) 
      END DO
  END SUBROUTINE dfactl
!***********************************************************************
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
      subroutine legend (plm,x,dfact,npt,lmax,m)
      IMPLICIT NONE
      INTEGER                          :: npt, lmax, m, maxfac
      REAL*8, DIMENSION(npt)           :: x
      REAL*8, DIMENSION(0:m)           :: dfact
      REAL*8, DIMENSION(npt,m:lmax)    :: plm
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
          plm(:,m) = dfact(m)
      ELSE
          fm=.5d+00*m
          plm(:,m) = dfact(m) * (1.d+00-x(:)*x(:))**fm
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
!***********************************************************************
!***********************************************************************
!deck sin_cos
!begin prologue     sin_cos
!date written       880721   (yymmdd)
!revision date      yymmdd   (yymmdd)
!keywords           sin_cos, link 2702, legendre functions
!author             schneider, barry (lanl)
!source             m2702
!purpose            sinmx and cosmx functions
!description        calculation of sin (m phi) and cos (m phi) functions
!references         
!                   The sin and cos functions for a given m value are computed on the    
!                   angular gridc    
!routines called
!end prologue       sin_cos
      subroutine sin_cos (sin_m,cos_m,x,npt,m)
      IMPLICIT NONE
      INTEGER                          :: npt, m
      REAL*8, DIMENSION(npt)           :: x
      REAL*8, DIMENSION(npt)           :: sin_m, cos_m
      REAL*8                           :: inv_sqrt_pi, inv_sqrt_2_pi
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      IF ( m == 0) THEN
           inv_sqrt_2_pi = one/sqrt(2.d0*pi)
           cos_m = inv_sqrt_2_pi
      ELSE
           inv_sqrt_pi = one/sqrt(pi)
           sin_m = inv_sqrt_pi*sin(m*x)
           cos_m = inv_sqrt_pi*cos(m*x)
      END IF
  END SUBROUTINE sin_cos
!***********************************************************************
!deck Two_Mat
!***begin prologue     Two_Mat
!***date written       030131   (yymmdd)
!***revision date               (yymmdd)
!***keywords           two electron radial integrals
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate two electron radial integrals
!***                   in a dvr representation by solving the
!***                   poisson equation for the inner integral
!***                   and then applying quadrature. 

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Two_Mat

  SUBROUTINE Two_Mat(v_two,ke,ang_pot,pt,wt,nphy)
  IMPLICIT NONE
  INTEGER                                :: nphy
  REAL*8,   DIMENSION(nphy)              :: pt, wt, ang_pot
  REAL*8,   DIMENSION(nphy,nphy)         :: ke
  REAL*8,   DIMENSION(nphy_tri,0:l_max)  :: v_two
  REAL*8,   DIMENSION(:,:), ALLOCATABLE  :: t
  REAL*8,   DIMENSION(:,:), ALLOCATABLE  :: work
  REAL*8                                 :: scale, ptlst 
  INTEGER,  DIMENSION(:), ALLOCATABLE    :: ipvt
  INTEGER                                :: l, lval, i, j, info, count
!
  scale=-2.d0*mass
  ALLOCATE(t(nphy-1,nphy-1),work(nphy,5),ipvt(nphy))
  DO l=0,l_max
!
!    put the second derivative sub-matrix into t and then scale it
!
     t(:,:)=ke(1:nphy-1,1:nphy-1)
     t = scale*t
     lval=l*(l+1)
!
!    add in the angular part
!
     DO i=1,nphy-1
        t(i,i) = t(i,i) - lval*ang_pot(i)
     END DO
!
!    invert t
!
     CALL dsytrf('l',nphy-1,t,nphy-1,ipvt,work,5*nphy,info)
     CALL dsytri('l',nphy-1,t,nphy-1,ipvt,work,info)
!
!    calculate the integrals by using the computed t to solve
!    the Poisson equation with homogeneous boundary conditions
!    and then add on the solution of the homogeneous Laplace
!    equation to get the result required.  The last point is special
!    and it contains no contribution from the sum.
!
     work(:,1)=pt**l
     work(:,2)=1.d0/(pt*sqrt(wt))
     lval=2*l+1
     ptlst=1.d0/(pt(nphy)**lval)
     count=0
     DO i=1,nphy-1
        DO j=1,i
           count=count+1
           v_two(count,l) = -lval*t(i,j)*work(i,2)*work(j,2) +   &
                             work(i,1)*work(j,1)*ptlst
        end DO
     END DO
     DO j=1,nphy
        count=count+1
        v_two(count,l) = work(nphy,1)*work(j,1)*ptlst
     END DO
  END DO
  DEALLOCATE(t,work,ipvt)
END SUBROUTINE Two_Mat
END MODULE CC_Prop_Module

