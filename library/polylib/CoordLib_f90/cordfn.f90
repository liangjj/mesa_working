! \usepackage{amsmath}
! \setkeys{Gin}{width=\linewidth}
! \title{Main Subroutine for Coordinate Functions}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck cordfn.f
!***begin prologue     cordfn
!***date written       000702   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr
!***
!***author             schneider, b. i.(nsf)
!***source             cordfn
!***purpose            get points and weights of Gauss quadratures for
!***                   general weight functions and then compute the
!***                   coordinate eigenfunctions.

!***references         see papers and notes appended.

!***routines called    iosys, util and mdutil
!***end prologue       cordfn
!
!     This is the main library routine to compute the orthogonal and
!     coordinate DVR functions for general weight functions.  The approach
!     is to use a reference quadrature to compute the $\alpha$ and $\beta$
!     recursion coefficients satisfied by the orthogonal polynonials.  The
!     three term recursion relationship is then diagonalized to obtain the
!     generalized points and weights.  The eigenvectors are used to transform
!     from the orthogonal polynomial to coordinate representation.  First
!     and second derivatives are also computed by taking the derivative of the
!     recursion relationship.
!
  SUBROUTINE cordfn(a,b,pt,wpt,cp,dcp,ddcp,wtt,dwtt,ddwtt,x,  &
                    alpha,beta,refalf,refbet,edge,ptfix,fix,  &
                    n,nq,npt,nfix,coord,stndrd,parity,typwt,  &
                    refwt,typint,prn)
  USE dvr_global,   ONLY   : inp, iout
  INTEGER                                :: n, nq, npt, nfix
  REAL*8, DIMENSION(n+1)                 :: a, b, pt, wpt 
  REAL*8, DIMENSION(npt,n+1)             :: cp, dcp, ddcp
  REAL*8, DIMENSION(npt)                 :: wtt, dwtt, ddwtt, x
  REAL*8                                 :: alpha, beta, refalf, refbet
  REAL*8, DIMENSION(2)                   :: edge
  INTEGER                                :: ptfix
  LOGICAL, DIMENSION(2)                  :: fix
  CHARACTER (LEN=*)                      :: coord, parity, typwt, refwt, typint
  LOGICAL                                :: stndrd, dollar
  LOGICAL, DIMENSION(10)                 :: prn
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=320)                    :: card
  REAL*8                                 :: mu, alf, bet, alfmax, dum
  REAL*8, DIMENSION(2)                   :: endpts
  REAL*8, DIMENSION(:),   ALLOCATABLE    :: r, rwt, wtfn, scr, arg
  REAL*8, DIMENSION(:,:), ALLOCATABLE    :: ply, eigv
  REAL*8                                 :: scaprd
  endpts(1)=edge(1)
  endpts(2)=edge(2)
!     The desire is to compute accurate values of the orthogonal polynomials
!     and co-ordinate functions for $i=0,1,2....n$.  So, we need to find
!     quadratures which are based on the zeros of the next higher polynomial.
  nord = n + 1
  IF(stndrd) THEN
!  
!     If the weight function is a one of the classical weight functions the
!     points and weights are known analytically and after computing them we
!     go directly to getting the coordinate functions.
     CALL gaussr(typwt,nord,alpha,beta,nfix,ptfix,b,pt,wpt,mu)
     mu=( edge(2)-edge(1) )*mu*.5D0
     WRITE(iout,*) ' mu = ',mu
     CALL chnvar(pt,wpt,-1.d0,1.d0,edge(1),edge(2),a,nord)
     IF(prn(4)) THEN
        title='final nodes'
        CALL prntfm(title,pt,nord,1,nord,1,iout)
        title='final weights'
        CALL prntfm(title,wpt,nord,1,nord,1,iout)
     END IF
     CALL cpoly(cp,dcp,ddcp,pt,a,n,npt,parity,prn(5))
  ELSE
!  
!    Compute the non-classical recursion coefficients using the Lanczos method
!    based on computing the needed integrals via a reference quadrature.
     alf=alpha
     bet=beta
     ALLOCATE(r(nq),rwt(nq),wtfn(nq),scr(nq),ply(nq,0:nord+1),arg(nq))
!  
!        Here we compute the reference nodes and weights.
!  
     CALL gaussq(refwt,nq,refalf,refbet,0,ptfix,scr, r,rwt)
  
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
!        Compute the reference weight function and the ratio of the
!        actual to the  reference weight.
!  
     CALL genrwt(wtfn,dum,dum,r,typwt,alf,bet,.false.,edge,nq)
     CALL genrwt(scr,dum,dum,r,refwt,refalf,refbet,.false.,edge,nq)
!     CALL vdiv(xwt(wtfn),xwt(wtfn),xwt(scr),nq)
     wtfn=wtfn/scr
     IF(prn(2)) THEN
        title='ratio weight factor'
        CALL prntfm(title,wtfn,nq,1,nq,1,iout)
     END IF
!  
!        Generate the recursion coefficients numerically.
!  
     ply(:,0)=1.d0
     arg=ply(:,0)*ply(:,0)*rwt*wtfn
     mu=sum(arg,1)
     WRITE(iout,7) mu
     IF(parity == 'none') THEN
        arg=r
     ELSE
        IF( dollar('$v0('//coord//')',card,title,inp) ) THEN
             angmom=intkey(card,'angular-momentum',0,' ')
        END IF
        arg = r*r
        ieo=angmom/2
        tst = angmom - 2*ieo
        IF(tst == 1) THEN
           ply(:,0)=r
        END IF
        endpts(1)=endpts(1)*endpts(1)
        endpts(2)=endpts(2)*endpts(2)
     END IF
!  
!        Note that lancz will return the recursion coefficients on
!        the actual interval.  This is consistent with what is done for
!        the known cases.  The starting function is specified or taken
!        to be unity.
!  
     CALL lancz(ply,arg,a,b,rwt,wtfn,scr,nq,nord)
!       Get rid of the unneeded memory
     DEALLOCATE(r,rwt,wtfn,scr,ply,arg)  
!       Modify them if using a Lobatto quadrature
     CALL modab(a,b,nfix,endpts,nord)
     IF(prn(3)) THEN
        title='lanczos a coefficients'
        CALL prntfm(title,a,nord,1,nord,1,iout)
        title='lanczos b coefficients'
        CALL prntfm(title,b,n,1,n,1,iout)
     END IF
!       Get the points and weights and then compute the coordinate
!       functions
     ALLOCATE(eigv(nord,nord),scr(npt))
     pt=a
     scr=b
!        Generate the non-classical points and weights and the
!        transformation matrix from the orthogonal polynomials
!        to the co-ordinate functions.
     CALL genq(pt,scr,wpt,eigv,fix,endpts,mu,nord)
     IF(parity /= 'none') THEN
        CALL xsq2x(pt,edge,nord)
     END IF
     title='eigenvectors'
     CALL prntrm(title,eigv,nord,nord,nord,nord,iout)
     IF(prn(4)) THEN
        title='final nodes'
        CALL prntfm(title,pt,nord,1,nord,1,iout)
        title='final weights'
       CALL prntfm(title,wpt,nord,1,nord,1,iout)
     END IF
!  
!        Generate the needed functions at all required points.
  
     CALL cpoly(cp,dcp,ddcp,pt,scr,n,npt,parity,prn(5))
     CALL genwtf(pt,scr,wtt,dwtt,ddwtt,alpha,beta, edge,npt,typwt,prn(5))
     DEALLOCATE(eigv,scr)
  END IF
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
END SUBROUTINE cordfn
