c \documentclass{article}
c \usepackage{graphicx}
c \usepackage{dcolumn}
c \usepackage{amsmath}
c \setkeys{Gin}{width=\linewidth}
c \title{Drvply:A Diluted Down Version of Orpoly}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck drvp.f
c***begin prologue     drvp
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             drvply
c***purpose            Points, weights and coordinate functions for generalized
c***                   Gauss quadratures.
c***description        Lanczos recursion using a reference weight function
c***                   is used to generate the points and weights of Gauss quadratures
c***                   for generalized weight functions.  The eigenvector matrix
c***                   of the tridiagonal matrix is used to compute the
c***                   coordinate functions, their first and second derivatives.
c
c***references         see papers and notes appended.
c
c***routines called    iosys, util and mdutil
c***end prologue       drvp

c     This is the main library routine to compute the orthogonal and
c     coordinate DVR functions for general weight functions.  The approach
c     is to use a reference quadrature to compute the $\alpha$ and $\beta$
c     recursion coefficients satisfied by the orthogonal polynonials.  The
c     three term recursion relationship is then diagonalized to obtain the
c     generalized points and weights.  The eigenvectors are used to transform
c     from the orthogonal polynomial to coordinate representation.  First
c     and second derivatives are also computed by taking the derivative of the
c     recursion relationship.

      subroutine drvp(q,wt,p,dp,ddp,edge,typwt,cordsys,n,nq,prn)
      implicit integer (a-z)
      character*(*) cordsys, typwt
      character*24 typarg
      character*80 title
      logical prn, fix
      real*8 q, wt, p, dp, ddp, edge, endpts, mu
      real*8 ab, xwt, dum, ptfix
      dimension q(*), wt(*), p(*), dp(*), ddp(*)
      dimension edge(2), ptfix(2), ngot(2), prn(*), fix(2)
      dimension endpts(2)
      common/io/inp, iout      
      data ptfix / -1.d0, 1.d0 /
      data fix / .true., .true. /
      data typarg / 'linear' /
      pointer (pab,ab(1))
      pointer (pxwt,xwt(1))
      endpts(1)=edge(1)
      endpts(2)=edge(2)
      nord = n + 1
      if(cordsys.eq.'cylindrical') then
         typarg='quadratic'
      endif 
      if(typwt.eq.'legendre') then
         a = 1
         b = a + nord
         need=wpadti(b+nord)
         call getmem(need,pab,ngot(1),'ab',0)

c        Get the recursion coefficients on [-1.0,1.0]

         call recur(typwt,nord,0.d0,0.d0,2,ptfix,ab(a),ab(b),mu)

c        Convert the coefficients to the interval [xl,xr]

         call cnvtab(ab(a),ab(b),edge,nord)
         mu=( edge(2)-edge(1) )*mu*.5d0
         if(prn(1)) then
            title='lanczos a coefficients'
            call prntfm(title,ab(a),nord,1,nord,1,iout)
            title='lanczos b coefficients'
            call prntfm(title,ab(b),n,1,n,1,iout)
            write(iout,*) ' mu = ',mu
         endif         
      else

c     Compute the non-classical recursion coefficients using the Lanczos method
c     based on computing the needed integrals via a reference quadrature.

         a=1
         b=a+nord
         r=b+nord
         rwt=r+nq
         wtfn=rwt+nq
         scr=wtfn+nq
         ply=scr+nq
         arg=ply+nq*(nord+1)
         need=wpadti(arg+nq)
         call getmem(need,pab,ngot(1),'ab',0)

c        Here we compute the reference nodes and weights.

         call gaussq('legendre',nq,0.d0,0.d0,0,ptfix,
     1                ab(scr),ab(r),ab(rwt))

c        Convert the Gauss points and weights on [-1.,1.] to [edge(1),edge(2)]

         call cnvtpt(ab(r),ab(rwt),edge,nq)      
         if(prn(1)) then
            title='reference nodes'
            call prntfm(title,ab(r),nq,1,nq,1,iout)
            title='reference weights'
            call prntfm(title,ab(rwt),nq,1,nq,1,iout)
         endif

c        Compute the weight function.

         call genrwt(ab(wtfn),dum,dum,ab(r),typwt,0.d0,0.d0,
     1              .false.,edge,nq)
         if(prn(2)) then
            title='ratio weight factor'
            call prntfm(title,ab(wtfn),nq,1,nq,1,iout)
         endif

c        Generate the recursion coefficients numerically.

         if(typarg.eq.'quadratic') then
            call vmul(ab(arg),ab(r),ab(r),nq)
            endpts(1)=endpts(1)*endpts(1)
            endpts(2)=endpts(2)*endpts(2)
         else  
            call copy(ab(r),ab(arg),nq)
         endif

c        Note that lancz will return the recursion coefficients on
c        the actual interval.  This is consistent with what is done for
c        the known cases.

         call lancz(ab(ply),ab(arg),ab(a),ab(b),ab(rwt),ab(wtfn),
     1              mu,ab(scr),nq,nord)
         call modab(ab(a),ab(b),2,endpts,nord)
         if(prn(3)) then
            title='lanczos a coefficients'
            call prntfm(title,ab(a),nord,1,nord,1,iout)
            title='lanczos b coefficients'
            call prntfm(title,ab(b),n,1,n,1,iout)
         endif
      endif

c     Get the points and weights and then compute the coordinate
c     functions

      xp=1
      xdp=xp+nord*nord
      xddp=xdp+nord*nord
      eigv=xddp+nord*nord
      scr=eigv+nord*nord
      wtt=scr+nord
      dwtt=wtt+nord
      ddwtt=dwtt+nord
      need=wpadti(ddwtt+nord)
      call getmem(need,pxwt,ngot(2),'xwt',0)
      call copy(ab(a),q,nord)       
      call copy(ab(b),xwt(scr),nord)       

c     Generate the non-classical points and weights and the
c     transformation matrix from the orthogonal polynomials
c     to the co-ordinate functions.

      call genq(q,xwt(scr),wt,xwt(eigv),fix,endpts,mu,nord)
      if(typarg.eq.'quadratic') then      
         call toquad(q,edge,nord)
      endif
      if(prn(4)) then
         title='final nodes'
         call prntfm(title,q,nord,1,nord,1,iout)
         title='final weights'
         call prntfm(title,wt,nord,1,nord,1,iout)
      endif

c     Generate the needed functions at all required points.

      call poly(xwt(xp),xwt(xdp),xwt(xddp),p,dp,ddp,q,xwt(scr),
     1          ab(a),ab(b),xwt(eigv),mu,n,nord,typarg,prn(5))

c     The DVR library assumes that the polynomials are $\delta$
c     functions at the quadrature points.  Convert to this normalization

      call todelta(p,dp,ddp,wt,nord)
      call getmem(-ngot(1),pab,idum,'ab',idum)
      call getmem(-ngot(2),pxwt,idum,'xwt',idum)
      return
 1    format(/,20x,'orthogonal polynomial basis function code')
 2    format(/,1x,'number of runs     = ',i3,
     1       /,1x,'interval type      = ',a16,
     1       /,1x,'weight function    = ',a16)
 3    format(/,1x,'left boundary          = ',e15.8,
     1       /,1x,'right boundary         = ',e15.8,
     2       /,1x,'number of fixed points = ',i1)
 4    format(/,1x,'alpha = ',e15.8,
     1       /,1x,'beta = ',e15.8) 
 5    format(/,1x,'initial rys parameter = ',e15.8,
     2       /,1x,'rys stepsize          = ',e15.8)
 6    format(/,1x,'polynomial n                 = ',i3,
     1       /,1x,'size of reference quadrature = ',i3,
     2       /,1x,'reference weight function    = ',a32,
     3       /,1x,'reference alpha              = ',e15.8,
     4       /,1x,'reference beta               = ',e15.8,
     5       /,1x,'rys alpha                    = ',e15.8)
      end
