! \documentclass{article}
! \usepackage{graphicx}
! \usepackage{dcolumn}
! \usepackage{amsmath}
! \setkeys{Gin}{width=\linewidth}
! \title{Drvply:A Diluted Down Version of Cordfn}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
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

  SUBROUTINE drvply(q,wt,p,dp,ddp,edge,typwt,angmom,nord,nq,reg_number)
  USE dvr_global,       ONLY  : inp, iout, l_val, m_val, nfix, nreg, fix
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nord
  REAL*8, DIMENSION(nord)                :: q, wt
  REAL*8, DIMENSION(nord,nord)           :: p, dp, ddp
  REAL*8, DIMENSION(2)                   :: edge
  CHARACTER (LEN=*)                      :: typwt
  INTEGER                                :: angmom
  INTEGER                                :: n, nq, reg_number, n_fixed
  CHARACTER (LEN=80)                     :: title
  REAL*8, DIMENSION(2)                   :: endpts
  REAL*8, DIMENSION(2)                   :: ptfix
  REAL*8                                 :: mu
  REAL*8, DIMENSION(:), ALLOCATABLE      :: a, b
  REAL*8, DIMENSION(:), ALLOCATABLE      :: r, rwt, wtfn, scr, arg, scrat
  REAL*8, DIMENSION(:,:), ALLOCATABLE    :: ply, eigv
  REAL*8                                 :: dum
  INTEGER                                :: err, tst, ieo
  DATA ptfix / -1.d0, 1.d0 /
!
  endpts=edge
  ALLOCATE(a(nord),b(nord),stat=err)
  if(err /= 0) then
     call lnkerr('allocation error')
  end if
  IF (nreg == 1) THEN
      n_fixed=0
      IF(nfix == 1) THEN
         n_fixed=1       
         ptfix(1) = -1.d0
         IF(fix(2)) THEN
            ptfix(1) = 1.d0
         END IF
      ELSE IF(nfix == 2) THEN
         n_fixed=2      
         ptfix(1) = -1.d0
         ptfix(2) = 1.d0
      END IF
  ELSE
      IF(reg_number == 1) THEN 
         n_fixed=1       
         ptfix(1)=1.d0
         IF(fix(1)) THEN
            n_fixed=2
            ptfix(1)=-1.d0
            ptfix(2)=1.d0
         END IF
      ELSE IF(reg_number == nreg) THEN 
         n_fixed=1   
         ptfix(1)=-1.d0    
         IF(fix(2)) THEN
            n_fixed=2
            ptfix(1)=-1.d0
            ptfix(2)=1.d0
         END IF
      ELSE
         n_fixed=2
         ptfix(1)=-1.d0
         ptfix(2)=1.d0
      END IF
  END IF
  IF(typwt == 'one'.or.         &
     typwt == 'chebyshev-1'.or. &
     typwt == 'chebyshev-2'.or. &
     typwt == 'jacobi'.or.      &
     typwt == 'hermite'.or.     &
     typwt == 'legendre'.or.    &
     typwt == 'laguerre')  THEN
!
!     If the weight function is a one of the classical weight functions the
!     points and weights are known analytically and after computing them we
!     go directly to getting the coordinate functions.
!
     CALL gaussq(typwt,nord,0.d0,0.d0,n_fixed,ptfix,b,q,wt)
     CALL cnvtpt(q,wt,edge,nord)
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
     CALL genrwt(wtfn,dum,dum,r,typwt,0.d0,0.d0,.false.,edge,nq)
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
    CALL modab(a,b,n_fixed,endpts,nord)
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
  CALL cpoly(p,dp,ddp,q,nord-1,nord,prn(5))
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
END SUBROUTINE drvply
