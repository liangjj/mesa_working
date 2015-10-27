c \documentclass{article}
c \usepackage{graphicx}
c \usepackage{dcolumn}
c \usepackage{amsmath}
c \setkeys{Gin}{width=\linewidth}
c \title{Main Subroutine for DVR Functions}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck orpoly.f 
c***begin prologue     orpoly
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             orpoly
c***purpose            get points and weights of Gauss quadratures for
c***                   general weight functions and then compute the
c***                   coordinate eigenfunctions.
c
c***references         see papers and notes appended.
c
c***routines called    iosys, util and mdutil
c***end prologue       orpoly

c     This is the main library routine to compute the orthogonal and
c     coordinate DVR functions for general weight functions.  The approach
c     is to use a reference quadrature to compute the $\alpha$ and $\beta$
c     recursion coefficients satisfied by the orthogonal polynonials.  The
c     three term recursion relationship is then diagonalized to obtain the
c     generalized points and weights.  The eigenvectors are used to transform
c     from the orthogonal polynomial to coordinate representation.  First
c     and second derivatives are also computed by taking the derivative of the
c     recursion relationship.

      subroutine orpoly(pntab,pntfun,x,alpha,beta,refalf,refbet,edge,
     1                  ptfix,fix,n,nq,npt,nfix,coord,stndrd,typarg,
     2                  typwt,refwt,typint,prn)
c
      implicit integer (a-z)
      character*(*) typwt, refwt, coord
      character*80 title
      character*(*) typint, typarg
      logical prn, fix, stndrd
      real*8 edge, alpha, beta, mu
      real*8 alf, bet, refalf, refbet, alfmax
      real*8 ab, xwt, x, dum, endpts
      real*8 y
#ifdef DECPOINTER
      integer*8 pntab, pntfun
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pntab, pntfun
#endif SGIPOINTER
      dimension x(npt), prn(10), edge(2), endpts(2), fix(2), ngot(3)
      common/io/inp, iout      

c     Pointer pntab will return location of recursion coeffients and
c     pntfun the needed points, wights, function and derivatives.

      pointer (pab,ab(1))
      pointer (pxwt,xwt(1))
      pointer (py,y(1))
      endpts(1)=edge(1)
      endpts(2)=edge(2)

c     The desire is to compute accurate values of the orthogonal polynomials
c     and co-ordinate functions for $i=0,1,2....n$.  So, we need to find
c     quadratures which are based on the zeros of the next higher polynomial.

      nord = n + 1
      
      if(stndrd) then

c     If the weight function is a one of the classical weight functions the
c     recursion coefficients are known analytically.

         a = 1
         b = a + nord
         need=wpadti(b+nord)
         call getmem(need,pab,ngot(1),'ab',0)

c        Get the recursion coefficients on [-1.0,1.0]

         call recur(typwt,nord,alpha,beta,nfix,ptfix,ab(a),ab(b),mu)

c        Convert the coefficients to the interval [edge(1),edge(2)]

         call cnvtab(ab(a),ab(b),edge,nord)
         mu=( edge(2)-edge(1) )*mu*.5d0
         if(prn(3)) then
            title='lanczos a coefficients'
            call prntfm(title,ab(a),nord,1,nord,1,iout)
            title='lanczos b coefficients'
            call prntfm(title,ab(b),n,1,n,1,iout)
            write(iout,*) ' mu = ',mu
         endif         
      else

c     Compute the non-classical recursion coefficients using the Lanczos method
c     based on computing the needed integrals via a reference quadrature.

         alf=alpha
         bet=beta
         a=1
         b=a+nord
         need=wpadti(b+nord)
         call getmem(need,pab,ngot(1),'ab',0)
         r=1
         rwt=r+nq
         wtfn=rwt+nq
         scr=wtfn+nq
         ply=scr+nq
         arg=ply+nq*(nord+1)
         need=wpadti(arg+nq)
         call getmem(need,pxwt,ngot(2),'xwt',0)

c        Here we compute the reference nodes and weights.

         call gaussq(refwt,nq,refalf,refbet,0,ptfix,xwt(scr),
     1               xwt(r),xwt(rwt))

c        Convert the Gauss points and weights on [-1.,1.] to [edge(1),edge(2)]

         call cnvtpt(xwt(r),xwt(rwt),edge,nq)      
         if(prn(1)) then
            title='reference nodes'
            call prntfm(title,xwt(r),nq,1,nq,1,iout)
            title='reference weights'
            call prntfm(title,xwt(rwt),nq,1,nq,1,iout)
         endif

c        Compute the reference weight function and the ratio of the
c        actual to the  reference weight.

         call genrwt(xwt(wtfn),dum,dum,xwt(r),typwt,alf,bet,
     1               .false.,edge,nq)
         call genrwt(xwt(scr),dum,dum,xwt(r),refwt,refalf,refbet,
     1              .false.,edge,nq)
         call vdiv(xwt(wtfn),xwt(wtfn),xwt(scr),nq)
         if(prn(2)) then
            title='ratio weight factor'
            call prntfm(title,xwt(wtfn),nq,1,nq,1,iout)
         endif

c        Generate the recursion coefficients numerically.

         if(typarg.eq.'quadratic') then
            call vmul(xwt(arg),xwt(r),xwt(r),nq)
            endpts(1)=endpts(1)*endpts(1)
            endpts(2)=endpts(2)*endpts(2)
         else  
            call copy(xwt(r),xwt(arg),nq)
         endif

c        Note that lancz will return the recursion coefficients on
c        the actual interval.  This is consistent with what is done for
c        the known cases.  The starting function is specified or taken
c        to be unity.

         call lancz(xwt(ply),xwt(arg),ab(a),ab(b),xwt(rwt),xwt(wtfn),
     1              mu,xwt(scr),nq,nord)
   
c        Get rid of the unneeded memory

         call getmem(-ngot(2),pxwt,idum,'xwt',idum)

c        Modify them if using a Lobatto quadrature

         call modab(ab(a),ab(b),nfix,endpts,nord)
         if(prn(3)) then
            title='lanczos a coefficients'
            call prntfm(title,ab(a),nord,1,nord,1,iout)
            title='lanczos b coefficients'
            call prntfm(title,ab(b),n,1,n,1,iout)
         endif
      endif

c     Get the points and weights and then compute the coordinate
c     functions

      pt=1
      wpt=pt+nord
      cp=wpt+nord
      dcp=cp+nord*npt
      ddcp=dcp+nord*npt
      wtt=ddcp+nord*npt
      dwtt=wtt+npt
      ddwtt=dwtt+npt
      kmat=ddwtt+npt
      need=wpadti(kmat+npt*npt)
      call getmem(need,pxwt,ngot(2),'xwt',0)
      p=1
      dp=p+nord*npt
      ddp=dp+nord*npt
      eigv=ddp+nord*npt
      scr=eigv+nord*nord
      need=wpadti(scr+npt)
      call getmem(need,py,ngot(3),'junk',0)
      call copy(ab(a),xwt(pt),nord)       
      call copy(ab(b),y(scr),nord)       

c     Generate the non-classical points and weights and the
c     transformation matrix from the orthogonal polynomials
c     to the co-ordinate functions.

      call genq(xwt(pt),y(scr),xwt(wpt),y(eigv),fix,endpts,mu,nord)
      title='eigenvectors'
      call prntrm(title,y(eigv),nord,nord,nord,nord,iout)
      if(typarg.eq.'quadratic') then      
         call toquad(xwt(pt),edge,nord)
      endif
      if(prn(4)) then
         title='final nodes'
         call prntfm(title,xwt(pt),nord,1,nord,1,iout)
         title='final weights'
         call prntfm(title,xwt(wpt),nord,1,nord,1,iout)
      endif
      call iosys('write integer "generalized quadrature size" to rwf',
     1            1,nord,0,' ')
      call iosys('write real "generalized alpha coefficients" to rwf',
     1            nord,ab(a),0,' ')
      call iosys('write real "generalized beta coefficients" to rwf',
     1            nord,ab(a),0,' ')
      call iosys('write real "generalized gauss points" to rwf',
     1            nord,xwt(pt),0,' ')
      call iosys('write real "generalized gauss weights" to rwf',
     1            nord,xwt(wpt),0,' ')

c     Generate the needed functions at all required points.

      call poly(y(p),y(dp),y(ddp),xwt(cp),xwt(dcp),xwt(ddcp),
     1          xwt(pt),y(scr),ab(a),ab(b),y(eigv),mu,n,
     2          npt,typarg,prn(5))
      call genwtf(xwt(pt),y(scr),xwt(wtt),xwt(dwtt),xwt(ddwtt),
     1            alpha,beta,edge,npt,typwt,prn(5))
      call getmem(-ngot(3),py,idum,'junk',idum)
      pntab=pab
      pntfun=pxwt
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
