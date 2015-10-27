!***********************************************************************
! FEDVR_Polynomials_Module
!**begin prologue     FEDVR_Polynomials_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Calculate the kinetic and potential energy matrix elements
!***                  in a FEDVR basis
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      FEDVR_Polynomials_Module
!***********************************************************************
!***********************************************************************
                           MODULE FEDVR_Polynomials_Module
                           USE FEDVR_Global
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
                            INTERFACE LaGrange_Polynomials                             
                       MODULE PROCEDURE Polynomials,                    &
                                        Odd_Polynomials,                &
                                        Fourier_Polynomials             
                            END INTERFACE LaGrange_Polynomials
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Polynomials.f
!***begin prologue     Polynomials
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            1. calculate piecewise lobatto dvr functions and
!***                      their one-body matrices
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Polynomials

  SUBROUTINE Polynomials(grid,reg_poly)
  IMPLICIT NONE
  TYPE(coordinates)                   :: grid
  TYPE(functions),  DIMENSION(:)      :: reg_poly
  REAL(idp)                           :: dum
  INTEGER                             :: i
!
!
  DO  i=1,nreg
      Write(iout,2) i, npt(i)
!
!        calculate the sector functions and their derivatives.
!
      ALLOCATE( grid%reg_pt_wt(i)%qr(npt(i)),                               &
                grid%reg_pt_wt(i)%wtr(npt(i)),                              &
                grid%reg_poly(i)%pr(npt(i),npt(i)),                         &
                grid%reg_poly(i)%dpr(npt(i),npt(i)),                        &
                grid%reg_poly(i)%ddpr(npt(i),npt(i)),                       &
                grid%reg_pt_wt(i)%inv_sqrt_wtr(npt(i)) )              
      CALL drvply( grid%reg_pt_wt(i)%qr,                                    &
                   grid%reg_pt_wt(i)%wtr,                                   &
                   grid%reg_poly(i)%pr,                                     &
                   grid%reg_poly(i)%dpr,                                    &
                   grid%reg_poly(i)%ddpr,                                   &
                   edge(i),                                                 &
                   typwt,                                                   &
                   npt(i),                                                  &
                   npt(i),                                                  &
                   i)
  END DO
!
!                  To compute most of what is required, it is not necessary
!                  to construct anything else. Since the bridge functions 
!                  span two elements, one can define them at the grid points
!                  but their derivatives are discontinuous across the
!                  sector boundaries. The matrix elements can be constructed 
!                  entirely from re-normalized sector quantities. 
!
  IF ( nreg == 1) THEN
!
!                  Only one region.  No endpoint corrections required.
!
       i = 1
       Call ReGrid ( grid%reg_pt_wt(i)%wtr,                                 &
                     grid%reg_pt_wt(i)%inv_sqrt_wtr,                        &
                     dum,                                                   &
                     dum,                                                   &
                     npt(i),                                                &
                     i)               
  ELSE
!
!                  First region.  Correction  at right endpoint needed from 
!                  first function in region 2.
       i = 1
       Call ReGrid ( grid%reg_pt_wt(i)%wtr,                                 &
                     grid%reg_pt_wt(i)%inv_sqrt_wtr,                        &
                     dum,                                                   &
                     grid%reg_pt_wt(i+1)%wtr(1),                            &
                     npt(i),                                                &
                     i)        
!       
       DO i = 2, nreg - 1
!
!                  General case.  Put result from the previous region into the
!                  the left region and correct the right endpoint.
!
          Call ReGrid ( grid%reg_pt_wt(i)%wtr,                              &
                        grid%reg_pt_wt(i)%inv_sqrt_wtr,                     &
                        grid%reg_pt_wt(i-1)%wtr(npt(i-1)),                  &
                        grid%reg_pt_wt(i+1)%wtr(1),                         &
                        npt(i),                                             &
                        i)               
       END DO
!
!                  Last region.  Correct the left end point.
!
       i = nreg
          Call ReGrid ( grid%reg_pt_wt(i)%wtr,                              &
                        grid%reg_pt_wt(i)%inv_sqrt_wtr,                     &
                        grid%reg_pt_wt(i-1)%wtr(npt(i-1)),                  &
                        dum,                                                &
                        npt(i),                                             &
                        i)               
  END IF
!
!                  Normalize the functions using the weights.
!
  DO i = 1, nreg 
     Call Re_Poly(                                                          &
                  grid%reg_poly(i)%pr,                                      &
                  grid%reg_poly(i)%dpr,                                     &
                  grid%reg_poly(i)%ddpr,                                    &
                  grid%reg_pt_wt(i)%inv_sqrt_wtr,                           &
                  npt(i) )
  END DO
1 FORMAT(/,10x,'Calculate the Regional Basis Functions for Coordinate = ',a4)
2 FORMAT(/,10x,'Region = ',i3,2x,'Number of Polynomials = ', i3)
END SUBROUTINE Polynomials
!***********************************************************************
!***********************************************************************
!deck Odd_Polynomials.f
!***begin prologue     Odd_Polynomials
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            1. calculate piecewise lobatto dvr functions and
!***                      their one-body matrices
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Odd_Polynomials

  SUBROUTINE Odd_Polynomials(grid,reg_poly_odd)
  IMPLICIT NONE
  TYPE(coordinates)                    :: grid
  TYPE(odd_functions),  DIMENSION(:)   :: reg_poly_odd
  REAL(idp)                            :: dum
  INTEGER                              :: i
  INTEGER                              :: j
!
!
  DO  i=1,nreg
      Write(iout,2) i, npt(i)
!
!        calculate the sector functions and their derivatives for the m odd grid.
!
      ALLOCATE( grid%reg_poly_odd(i)%pr(npt(i),npt(i)),                      &
                grid%reg_poly_odd(i)%dpr(npt(i),npt(i)),                     &
                grid%reg_poly_odd(i)%ddpr(npt(i),npt(i)) )
      DO j = 1, npt(i)
         grid%reg_poly_odd(i)%pr(:,j)      =                                 &
         grid%reg_pt_wt(i)%inv_qr_fac(j) * grid%reg_poly(i)%pr(:,j) 
         grid%reg_poly_odd(i)%dpr(:,j)     =                                 &
         grid%reg_pt_wt(i)%inv_qr_fac(j) * grid%reg_poly(i)%dpr(:,j) 
         grid%reg_poly_odd(i)%ddpr(:,j)    =                                 &
         grid%reg_pt_wt(i)%inv_qr_fac(j) * grid%reg_poly(i)%ddpr(:,j) 
      END DO
  END DO
1 FORMAT(/,10x,'Calculate the Regional Basis Functions for Coordinate = ',a4)
2 FORMAT(/,10x,'Region = ',i3,2x,'Number of Polynomials = ', i3)
END SUBROUTINE Odd_Polynomials
!***********************************************************************
!***********************************************************************
!***deck Fourier_Polynomials.f                                                                          
!***begin prologue     Fourier_Polynomials.f
!***date written       000702   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr
!***
!***author             schneider, b. i.(nsf)
!***source             Fourier_Points_Weights
!***purpose            Points, and weights for fourier dvr
!***description
!***
!***
!***
!*** 
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Fourier_Polynomials
  SUBROUTINE Fourier_Polynomials(grid,reg_poly_fourier)
  IMPLICIT NONE
  TYPE(coordinates)                      :: grid
  TYPE(fourier_functions),  DIMENSION(:) :: reg_poly_fourier
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: k
  REAL(idp)                              :: fac
  CHARACTER (LEN=80)                     :: title
  ALLOCATE(grid%reg_poly_fourier(1))
  ALLOCATE(grid%reg_pt_wt(1)%qr(npt(1)),           &
           grid%reg_pt_wt(1)%wtr(npt(1)),          &
           grid%reg_poly_fourier(1)%pr(npt(1),npt(1)))
  box=edge(2)-edge(1)
  deltax=box/npt(1)
  grid%reg_pt_wt(1)%wtr(:) = deltax
  j = ( npt(1) - 1)/2
  fac = edge(1) + .5d0*box
  k=-j
  DO i=1,npt(1)
     grid%reg_pt_wt(1)%qr(i) = deltax * k + fac
     k=k+1
  END DO
  IF(prn(3)) THEN
     title='fourier points'
     CALL prntrm(title,grid%reg_pt_wt(1)%qr,n,1,n,1,iout)
     title='fourier weights'
     CALL prntrm(title,grid%reg_pt_wt(1)%wtr,n,1,n,1,iout)
  END IF
  grid%reg_poly_fourier(1)%pr(:,:) = zero
  fac=1.d0/sqrt(deltax)
  DO i=1,npt(1)
     grid%reg_poly_fourier(1)%pr(i,i) = fac
  END DO
  END SUBROUTINE Fourier_Polynomials
!***********************************************************************
!***********************************************************************
!deck ReGrid.f
!***begin prologue     ReGrid
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Modify the weights at the ends of the interval to
!***                   reflect that there are bridge functions present. 
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       ReGrid

  SUBROUTINE ReGrid(wt,inv_sqrt_wt,wt_left_end,wt_right_end,n,region)
  IMPLICIT NONE
  REAL(idp), DIMENSION (:)             :: wt
  REAL(idp), DIMENSION (:)             :: inv_sqrt_wt
  REAL(idp)                            :: wt_left_end
  REAL(idp)                            :: wt_right_end
  INTEGER                              :: n
  INTEGER                              :: region
!
!
  IF ( nreg == 1) THEN
       inv_sqrt_wt(1:n) = Sqrt ( 1.d0 / wt(1:n) )
       Return
  END IF
  IF ( region == 1) THEN
!
!      Modify the last weight and then get the inverse square roots.
!
       wt(n) = wt(n) + wt_right_end
       inv_sqrt_wt(1:n) = 1.d0 / sqrt ( wt(1:n) )     
  ELSE IF ( region == nreg ) THEN
!
!      
!
       wt(1) = wt_left_end
       inv_sqrt_wt(1:n) = 1.d0 / sqrt ( wt(1:n) )
  ELSE
!      
!      Modify the last weight and then get the inverse square roots.
!
       wt(1)  = wt_left_end
       wt(n)  = wt(n) + wt_right_end
       inv_sqrt_wt(1:n)  =  1.d0 / sqrt ( wt(1:n) )
  END IF
END SUBROUTINE ReGrid
!***********************************************************************
!***********************************************************************
!deck Re_Poly.f
!***begin prologue     Re_Poly
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            1. calculate piecewise lobatto dvr functions and
!***                      their one-body matrices
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Re_Poly

  SUBROUTINE Re_Poly(p,dp,ddp,inv_sqrt_wt,n)
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)          :: p
  REAL*8, DIMENSION(:,:)          :: dp
  REAL*8, DIMENSION(:,:)          :: ddp
  REAL*8, DIMENSION (:)           :: inv_sqrt_wt
  INTEGER                         :: n
  INTEGER                         :: i
!
!
  DO i = 1, n
       p(:, i ) =    p(:, i )  * inv_sqrt_wt(i)
      dp(:, i ) =   dp(:, i )  * inv_sqrt_wt(i)
     ddp(:, i ) =  ddp(:, i )  * inv_sqrt_wt(i)
  END DO
END SUBROUTINE Re_Poly
!***********************************************************************
!***********************************************************************
!deck drvply.f
!***begin prologue     drvply
!***date written       000702   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr
!***
!***author             schneider, b. i.(nsf)
!***source             drvply
!***purpose            Points, weights and coordinate functions for generalized
!***                   Gauss quadratures.
!***description        Lanczos recursion using a reference weight function
!***                   is used to generate the points and weights of Gauss quadratures
!***                   for generalized weight functions.  The eigenvector matrix
!***                   of the tridiagonal matrix is used to compute the
!***                   coordinate functions, their first and second derivatives.

!***references         see papers and notes appended.

!***routines called    iosys, util and mdutil
!***end prologue       drvply

!     This is the main library routine to compute the orthogonal and
!     coordinate DVR functions for general weight functions.  The approach
!     is to use a reference quadrature to compute the $\alpha$ and $\beta$
!     recursion coefficients satisfied by the orthogonal polynonials.  The
!     three term recursion relationship is then diagonalized to obtain the
!     generalized points and weights.  The eigenvectors are used to transform
!     from the orthogonal polynomial to coordinate representation.  First
!     and second derivatives are also computed by taking the derivative of the
!     recursion relationship.

 SUBROUTINE drvply(q,wt,p,dp,ddp,edge,typwt,nord,nq,reg_number)
  IMPLICIT NONE
  INTEGER                                :: nord
  REAL(idp), DIMENSION(:)                :: q
  REAL(idp), DIMENSION(:)                :: wt
  REAL(idp), DIMENSION(:,:)              :: p
  REAL(idp), DIMENSION(:,:)              :: dp
  REAL(idp), DIMENSION(:,:)              :: ddp
  REAL(idp), DIMENSION(2)                :: edge
  CHARACTER (LEN=*)                      :: typwt
  INTEGER                                :: n
  INTEGER                                :: nq
  INTEGER                                :: reg_number
  INTEGER                                :: n_fixed
  CHARACTER (LEN=80)                     :: title
  REAL(idp), DIMENSION(2)                :: endpts
  REAL(idp), DIMENSION(2)                :: ptse
  REAL(idp), DIMENSION(2)                :: ptfix
  REAL(idp)                              :: mu
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: a
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: b
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: r
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: rwt
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: wtfn
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: scr
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: arg
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: scrat
  REAL(idp), DIMENSION(:,:), ALLOCATABLE :: ply
  REAL(idp), DIMENSION(:,:), ALLOCATABLE :: eigv
  REAL(idp)                              :: dum
  INTEGER                                :: err
  INTEGER                                :: tst
  INTEGER                                :: ieo
  DATA ptfix / -1.d0, 1.d0 /
!
  ALLOCATE(a(nord),b(nord),stat=err)
  if(err /= 0) then
     call lnkerr('allocation error')
  end if
!
! If the weight function is a one of the classical weight functions the
! points and weights are known analytically and after computing them we
!  go directly to getting the coordinate functions.
!
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
! For finite regions
!
  write(iout,8) typwt
  IF(typwt == 'one'.or.         &
     typwt == 'chebyshev-1'.or. &
     typwt == 'chebyshev-2'.or. &
     typwt == 'jacobi'.or.      &
     typwt == 'spherical'.or.   &
     typwt == 'cylindrical'.or. &
     typwt == 'legendre') THEN
     CALL gaussq(typwt,nord,0.d0,0.d0,n_fixed,ptfix,b,q,wt)
     CALL cnvtpt(q,wt,edge,nord)
!
!  For Infinite Regions
!
  ELSE IF( typwt == 'hermite') THEN
     n_fixed = 0
     CALL gaussq(typwt,nord,0.d0,0.d0,n_fixed,ptfix,b,q,wt)
  ELSE IF(typwt == 'laguerre')  THEN
     n_fixed=0
     IF(nfix == 1) THEN
        n_fixed=1       
        ptfix(1) = 1.0d-20
        ptse(1) = ptfix(1)    
        CALL gaussq(typwt,nord,0.d0,0.d0,n_fixed,ptfix,b,q,wt)
     END IF
  ELSE
!  
!     Compute the non-classical recursion coefficients using the Lanczos method
!     based on computing the needed integrals via a reference quadrature.
!  
     ALLOCATE(r(nq),rwt(nq),wtfn(nq),scr(nq),ply(nq,0:nord),arg(nq),stat=err)
     if(err /= 0) then
        call lnkerr('allocation error')
     end if
!  
!        Here we compute the reference nodes and weights.
!  
     CALL gaussq('legendre',nq,0.d0,0.d0,0,ptfix,scr,r,rwt)
!  
!        Convert the Gauss points and weights on [-1.,1.] to [edge(1),edge(2)]
!
     CALL cnvtpt(r,rwt,edge,nq)
     IF(prn(1)) THEN
        title='reference nodes'
        CALL prntfm(title,r,nq,1,nq,1,iout)
        title='reference weights'
        CALL prntfm(title,rwt,nq,1,nq,1,iout)
     END IF
!  
!        Compute the weight function.
     CALL genrwt(wtfn,r,typwt,0.d0,0.d0,.false.,nq)
     IF(prn(2)) THEN
        title='ratio weight factor'
        CALL prntfm(title,wtfn,nq,1,nq,1,iout)
     END IF
!  
!        Generate the recursion coefficients numerically.
!        Initialize the first function and normalize.
!  
     ply(:,0) = 1.d0
     arg=ply(:,0)*ply(:,0)*rwt*wtfn
     mu=sum(arg,1)
     WRITE(iout,7) mu
!  
!        Note that lancz will return the recursion coefficients on
!        the actual interval.  This is consistent with what is done for
!        the known cases.
!  
    arg = r
    CALL lancz(ply,arg,a,b,rwt,wtfn,scr,nq,nord)
    CALL modab(a,b,n_fixed,ptse,nord)
    IF(prn(3)) THEN
        title='lanczos a coefficients'
        CALL prntfm(title,a,nord,1,nord,1,iout)
        title='lanczos b coefficients'
        CALL prntfm(title,b,nord-1,1,nord-1,1,iout)
    END IF
!  
!        Get the points and weights and then compute the coordinate
!        functions
!  
    ALLOCATE(eigv(nord,nord),scrat(nord),stat=err)
    if(err /= 0) then
       call lnkerr('allocation error')
    end if
    q = a
    scrat=b
!  
!        Generate the non-classical points and weights and the
!        transformation matrix from the orthogonal polynomials
!        to the co-ordinate functions.
!
    CALL genq(q,scrat,wt,eigv,fix,endpts,mu,nord)
!    IF(parity /= 'none') THEN
!       CALL xsq2x(q,edge,nord)
!    END IF
    DEALLOCATE(r,rwt,wtfn,scr,ply,arg,eigv,scrat)
  END IF
  IF(prn(4)) THEN
     title='final nodes'
     CALL prntfm(title,q,nord,1,nord,1,iout)
     title='final weights'
     CALL prntfm(title,wt,nord,1,nord,1,iout)
  END IF
!  
! Generate the needed functions at all required points.
!  
!  CALL cpoly(p,dp,ddp,q,a,nord-1,nord,parity,prn(5))
  CALL cpoly(p,dp,ddp,q,a,nord-1,nord,prn(5))
!  
! The DVR library assumes that the polynomials are $\delta$
! functions at the quadrature points.  Convert to this normalization
!
  DEALLOCATE(a,b)  
!
1    FORMAT(/,20X,'orthogonal polynomial basis function code')
2    FORMAT(/,1X,'number of runs     = ',i3,  &
    /,1X,'interval type      = ',a16, /,1X,'weight function    = ',a16)
3    FORMAT(/,1X,'left boundary          = ',e15.8,  &
    /,1X,'right boundary         = ',e15.8, /,1X,'number of fixed points = ',i1)
4    FORMAT(/,1X,'alpha = ',e15.8, /,1X,'beta = ',e15.8)
5    FORMAT(/,1X,'initial rys parameter = ',e15.8,  &
    /,1X,'rys stepsize          = ',e15.8)
6    FORMAT(/,1X,'polynomial n                 = ',i3,  &
    /,1X,'size of reference quadrature = ',i3,  &
    /,1X,'reference weight function    = ',a32,  &
    /,1X,'reference alpha              = ',e15.8,  &
    /,1X,'reference beta               = ',e15.8,  &
    /,1X,'rys alpha                    = ',e15.8)
7    FORMAT(/,1X,'weight integral = ',e15.8)
8    FORMAT(/,1X,'Generating Point and Weights for a Weight Function = ',a16)
END SUBROUTINE drvply
!***********************************************************************
!***********************************************************************
!deck cnvtpt.f
  SUBROUTINE cnvtpt(pt,wt,endpts,n)
  USE input_output
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
!***********************************************************************
!***********************************************************************
!deck genq.f
  SUBROUTINE genq( a, b, w, z, fix, edge, muzero, n)
  IMPLICIT NONE
!        a, b     on input the recursion coefficients; will be destroyed
!     output parameters (both arrays of length n)
!        a        the desired points
!        w        the desired weights
  INTEGER                                :: n
  REAL(idp), DIMENSION(n)                :: a
  REAL(idp), DIMENSION(n)                :: b
  REAL(idp), DIMENSION(n)                :: w
  REAL(idp), DIMENSION(n,n)              :: z
  LOGICAL, DIMENSION(2)                  :: fix
  REAL(idp), DIMENSION(2)                :: edge
  REAL(idp)                              :: muzero
  REAL(idp)                              :: sumwt
  INTEGER                                :: i
  INTEGER                                :: ierr
!     the method used is a ql-type method with origin shifting
  z(:,1)=b
  b(2:n)=z(1:n-1,1)
  w=0.d0
  z=0.d0
  DO  i=1,n
      z(i,i)=1.d0
  END DO
  CALL imtql2 (n, n, a, b, z, ierr)
  w=muzero*z(1,:)*z(1,:)
  sumwt=sum(w,1)
  IF(fix(1)) THEN
     a(1)=edge(1)
  END IF
  IF(fix(2)) THEN
     a(n)=edge(2)
  END IF
  WRITE(iout,1) sumwt
1    FORMAT(/,1X,'sum of the weights = ',e15.8)
END SUBROUTINE genq
!***********************************************************************
!***********************************************************************
!deck genrwt.f
!***begin prologue     genrwt
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           reference weight
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate weight functions and derivatives.
!***
!***description
!***
!***
!***references
!***routines called
!***end prologue       genrwt
  SUBROUTINE genrwt(rwt,pt,wt,alpha,beta,deriv,n,drwt,ddrwt)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL(idp), DIMENSION(:)                :: rwt
  REAL(idp), DIMENSION(:)                :: pt
  CHARACTER (LEN=*)                      :: wt
  REAL(idp)                              :: alpha
  REAL(idp)                              :: beta
  LOGICAL                                :: deriv
  REAL(idp), OPTIONAL, DIMENSION(:)      :: drwt
  REAL(idp), OPTIONAL, DIMENSION(:)      :: ddrwt
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: fac1
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: fac2
!     Weight functions and their derivatives are computed for a variety of
!     cases
  IF(wt == 'one'.or.       &
     wt == 'legendre'.or.  &
     wt == 'spherical'.or. &
     wt == 'cylindrical') THEN
     rwt = 1.d0
     IF(deriv) THEN
        drwt=0.d0
        ddrwt=0.d0
     END IF
  ELSE IF(wt == 'r') THEN
     rwt=pt
     IF(deriv) THEN
        drwt=1.d0
        ddrwt=0.d0
     END IF
  ELSE IF(wt == 'rr') THEN
     rwt=pt*pt 
     IF(deriv) THEN
       drwt=2.d0*pt
       ddrwt=2.d0
     END IF
  ELSE IF(wt == 'hermite') THEN
    rwt=EXP(-pt*pt)
    IF(deriv) THEN
       drwt = -2.d0*rwt*pt
       ddrwt= ( -2.d0 + 4.d0*pt*pt ) * rwt
    END IF
  ELSE IF(wt == 'chebyshev-1') THEN
    rwt= SQRT ( 1.d0/( 1.d0 - pt*pt ) )
    IF(deriv) THEN
       ddrwt = rwt * rwt * rwt
       drwt = pt * ddrwt
       ddrwt = ddrwt +  3.d0*pt*pt*ddrwt*rwt*rwt
    END IF
  ELSE IF(wt == 'chebyshev-2') THEN
    rwt= SQRT (  1.d0 -pt*pt )
    IF(deriv) THEN
       drwt = - pt/rwt
       ddrwt = ( - 1.d0 + pt*drwt/rwt )/rwt
    END IF
  ELSE IF(wt == 'laguerre') THEN
    rwt= pt**alpha*EXP(-pt)
    IF(deriv) THEN
       ALLOCATE(fac1(n),fac2(n))
       fac1 =  - pt + alpha/pt 
       fac2 =  1.d0 + alpha/(pt*pt) 
       drwt = fac1 * rwt
       ddrwt = fac1 * drwt - fac2 * rwt
       DEALLOCATE(fac1,fac2)
    END IF
  ELSE IF(wt == 'jacobi') THEN
    rwt= (1.d0-pt)**alpha * (1.d0+pt)**beta
    IF(deriv) THEN
       ALLOCATE(fac1(n),fac2(n))
       fac1= -alpha/(1.d0-pt) 
       fac2=beta/(1.d0+pt) 
       drwt = (fac1 + fac2 )*rwt
       ddrwt = ( fac1 + fac2 )*drwt
       fac1=fac1/(1.d0-pt)
       fac2=fac2/(1.d0+pt)
       ddrwt = ddrwt + ( fac1 + fac2 )*rwt
       DEALLOCATE(fac1,fac2)
    END IF
  ELSE IF(wt == 'rys') THEN
    rwt=EXP(-alpha*pt*pt)
    IF(deriv) THEN
      ALLOCATE(fac1(n))
      fac1 = -2.d0*pt*alpha
      drwt = fac1*rwt
      ddrwt = fac1*drwt - 2.d0*alpha*rwt
      DEALLOCATE(fac1)
    END IF
  ELSE
    rwt=1.d0
  END IF
END SUBROUTINE genrwt
!***********************************************************************
!***********************************************************************
!deck lancz.f
!***begin prologue     lancz
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           eigenvalues, eigenvectors
!***author             schneider, barry (nsf)
!***source
!***purpose            generate orthogonal polynomials using recursion.
!***                   using the generated alpha and beta coefficients find
!***                   the points and weights of the generalized gauss
!***                   quadrature.
!***
!***description
!***references
!***routines called
!***end prologue       lancz
!\begin{eqnarray}
!  \beta_{j} P_{j}(x) &=& \big ( x - \alpha_{j} \big ) P_{j-1}(x) - \beta_{j-1} P_{j-2}(x)
!                             \\ \nonumber
!  \alpha_{j} &=& \langle P_{j-1} \mid x \mid P_{j-1} \rangle \\ \nonumber
!  \beta_{j} &=& \langle P_{j-1} \mid x \mid P_{j} \rangle \\ \nonumber
!  \beta_{j} P_{j}^{\prime}(x) &=& \big ( x - \alpha_{j} \big ) P_{j-1}^{\prime}(x)
!                             - \beta_{j-1} P_{j-2}^{\prime}(x) + P_{j-1}(x)
!                                 \\ \nonumber
!  \beta_{j} P_{j}^{\prime \prime}(x) &=& \big ( x - \alpha_{j} \big )
!                                      P_{j-1}^{\prime \prime}(x)
!                            - \beta_{j-1} P_{j-2}^{\prime \prime}(x)
!                            + 2 P_{j-1}^{\prime}(x)
!\end{eqnarray}
  SUBROUTINE lancz(v,arg,a,b,wt,refwt,scr,n,iter)
  IMPLICIT NONE
  INTEGER                                :: n
  INTEGER                                :: iter
  REAL(idp), DIMENSION(n,0:iter)         :: v
  REAL(idp), DIMENSION(n)                :: arg
  REAL(idp), DIMENSION(n)                :: wt
  REAL(idp), DIMENSION(n)                :: refwt
  REAL(idp), DIMENSION(n)                :: scr
  REAL(idp), DIMENSION(iter)             :: a
  REAL(idp), DIMENSION(iter)             :: b
  REAL(idp)                              :: anorm
  INTEGER                                :: i
  REAL(idp), DIMENSION(:), ALLOCATABLE   :: vtmp
  ALLOCATE(vtmp(n))
!     first vector is input.  normalize.
  vtmp=v(:,0)*v(:,0)*wt*refwt
  anorm=SQRT(1.d0/sum(vtmp,1))
  v(:,0)=anorm*v(:,0)
!     we now have the first function
  IF (iter > 1) THEN
!         form argument times the first function
      scr=arg*v(:,0)
!         calculate a(1)
      vtmp=v(:,0)*scr*wt*refwt
      a(1)=sum(vtmp,1)
!         form mat times the first function - a(1) times the first function
!         and store it in the next polynomial
      v(:,1)=scr - a(1)*v(:,0)
!         calculate b(1)
      vtmp=v(:,1)*v(:,1)*wt*refwt
      b(1)=SQRT( sum(vtmp,1) )
!         normalize the second polynomial
      v(:,1)=(1.d0/b(1))*v(:,1)
  END IF
  IF (iter > 2) THEN
      DO  i=2,iter
!            multiply the last calculated polynomial by mat
          scr=arg*v(:,i-1)
!            orthogonalize it to the two previous polynomials
!            calculating a(i) as we go
          vtmp=v(:,i-1)*scr*wt*refwt
          a(i)=sum(vtmp,1)
          v(:,i)=scr - a(i)*v(:,i-1) - b(i-1)*v(:,i-2)
!            calculate b(i)
          vtmp=v(:,i)*v(:,i)*wt*refwt       
          b(i) = SQRT(sum(vtmp,1))
!            normalize the polynomial and we are done
          v(:,i)=(1.d0/b(i))*v(:,i)
      END DO
  END IF
  DEALLOCATE(vtmp)
END SUBROUTINE lancz
!***********************************************************************
!**********************************************************************
!deck modab.f
  SUBROUTINE modab(a,b,kpts,endpts,n)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL(idp), DIMENSION(n)                :: a
  REAL(idp), DIMENSION(n)                :: b
  INTEGER                                :: kpts
  REAL(idp), DIMENSION(2)                :: endpts
  REAL(idp)                              :: gbslve
  REAL(idp)                              :: gam
  REAL(idp)                              :: t1
  IF (kpts == 0) THEN
      RETURN
  ELSE IF (kpts == 1) THEN
!         only a(n) must be changed
      a(n) =gbslve(endpts(1), n, a, b)*b(n-1)**2 + endpts(1)
      RETURN
  ELSE IF(kpts == 2) THEN
!         a(n) and b(n-1) must be recomputed
      gam =gbslve(endpts(1), n, a, b)
      t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, a, b) - gam))
      b(n-1) =  SQRT(t1)
      a(n) = endpts(1) + gam*t1
      RETURN
  END IF
END SUBROUTINE modab
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

  SUBROUTINE cpoly(cp,dcp,ddcp,pt,arg,n,npt,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                             :: n, npt
  REAL(idp), DIMENSION(npt,0:n)          :: cp, dcp, ddcp
  REAL(idp), DIMENSION(npt)              :: pt, arg
  LOGICAL                             :: prn
  REAL(idp)                              :: fac, fac2, tmp
  CHARACTER (LEN=80)                  :: title
  arg = pt
  CALL lgngr(cp,dcp,ddcp,pt,pt,npt,npt,.false.,'all')
  IF(prn) THEN
     title='coordinate function'
     CALL prntfm(title,cp(1,0),npt,n+1,npt,n+1,iout)
     title='first derivative of coordinate function'
     CALL prntfm(title,dcp(1,0),npt,n+1,npt,n+1,iout)
     title='second derivative of coordinate function'
     CALL prntfm(title,ddcp(1,0),npt,n+1,npt,n+1,iout)
   END IF
END SUBROUTINE cpoly
!***********************************************************************
!***********************************************************************
           END MODULE FEDVR_Polynomials_Module
!***********************************************************************
!***********************************************************************
