MODULE TESTMOD
                            CONTAINS
 SUBROUTINE drvply(q,wt,p,dp,ddp,edge,typwt,nord,nq,reg_number)
  USE dvr_global,       ONLY  : inp, iout, l_val, m_val, nfix, nreg, fix
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nord
  REAL*8, DIMENSION(:)                :: q, wt
  REAL*8, DIMENSION(:,:)           :: p, dp, ddp
  REAL*8, DIMENSION(2)                   :: edge
  CHARACTER (LEN=*)                      :: typwt
!  CHARACTER (LEN=*)                      :: parity
!  INTEGER                                :: angmom
  INTEGER                                :: n, nq, reg_number, n_fixed
  CHARACTER (LEN=80)                     :: title
  REAL*8, DIMENSION(2)                   :: endpts
  REAL*8, DIMENSION(2)                   :: ptfix
  REAL*8                                 :: mu
  REAL*8, DIMENSION(:), ALLOCATABLE      :: a, b
  REAL*8, DIMENSION(:), ALLOCATABLE      :: r, rwt, wtfn, scr, arg, scrat
  REAL*8, DIMENSION(:,:), ALLOCATABLE    :: ply, eigv
  REAL*8, DIMENSION(1)                                 :: dum
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
END SUBROUTINE drvply
! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Reference Weight Function}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
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
  SUBROUTINE genrwt(rwt,drwt,ddrwt,pt,wt,alpha,beta,deriv,edge,n)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(:)                   :: rwt, drwt, ddrwt, pt
  CHARACTER (LEN=*)                      :: wt
  REAL*8                                 :: alpha, beta
  LOGICAL                                :: deriv
  REAL*8, DIMENSION(2)                   :: edge
  REAL*8, DIMENSION(:), ALLOCATABLE      :: fac1, fac2
!     Weight functions and their derivatives are computed for a variety of
!     cases
  IF(wt == 'one'.or.wt == 'legendre') THEN
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
       ddrwt = rwt*rwt*rwt
       drwt = pt*ddrwt
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
END MODULE TESTMOD
