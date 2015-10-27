c \documentclass{article}
c \usepackage{graphicx}
c \usepackage{dcolumn}
c \usepackage{amsmath}
c \setkeys{Gin}{width=\linewidth}
c \title{Main Subroutine for Coordinate Functions}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck cordfn.f 
c***begin prologue     cordfn
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             cordfn
c***purpose            get points and weights of Gauss quadratures for
c***                   general weight functions and then compute the
c***                   coordinate eigenfunctions.
c
c***references         see papers and notes appended.
c
c***routines called    iosys, util and mdutil
c***end prologue       cordfn

c     This is the main library routine to compute the orthogonal and
c     coordinate DVR functions for general weight functions.  The approach
c     is to use a reference quadrature to compute the $\alpha$ and $\beta$
c     recursion coefficients satisfied by the orthogonal polynonials.  The
c     three term recursion relationship is then diagonalized to obtain the
c     generalized points and weights.  The eigenvectors are used to transform
c     from the orthogonal polynomial to coordinate representation.  First
c     and second derivatives are also computed by taking the derivative of the
c     recursion relationship.

      subroutine cordfn(a,b,pt,wpt,cp,dcp,ddcp,wtt,dwtt,ddwtt,x,
     1                  alpha,beta,refalf,refbet,edge,ptfix,fix,
     2                  n,nq,npt,nfix,coord,stndrd,parity,typwt,
     3                  refwt,typint,prn)
c
      implicit integer (a-z)
      character*(*) typwt, refwt, coord
      character*80 title
      character*320 card
      character*(*) typint, parity
      logical prn, fix, stndrd, dollar
      real*8 a, b, pt, wpt, cp, dcp, ddcp, wtt, dwtt, ddwtt
      real*8 edge, alpha, beta, mu
      real*8 alf, bet, refalf, refbet, alfmax
      real*8 x, dum, endpts
      real*8 xwt, y
      real*8 scaprd
      dimension a(n+1), b(n+1), pt(n+1), wpt(n+1)
      dimension cp(npt,n+1), dcp(npt,n+1), ddcp(npt,n+1)
      dimension wtt(npt), dwtt(npt), ddwtt(npt)
      dimension x(npt), prn(10), edge(2), endpts(2), fix(2), ngot(2)
      common/io/inp, iout      
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
c     points and weights are known analytically and after computing them we
c     go directly to getting the coordinate functions.    

         call gaussr(typwt,nord,alpha,beta,nfix,ptfix,b,pt,wpt,mu)
         mu=( edge(2)-edge(1) )*mu*.5d0
         write(iout,*) ' mu = ',mu
         call chnvar(pt,wpt,-1.d0,1.d0,edge(1),edge(2),a,nord)
         if(prn(4)) then
            title='final nodes'
            call prntfm(title,pt,nord,1,nord,1,iout)
            title='final weights'
            call prntfm(title,wpt,nord,1,nord,1,iout)
         endif
         call cpoly(cp,dcp,ddcp,pt,a,n,npt,parity,prn(5))
      else

c     Compute the non-classical recursion coefficients using the Lanczos method
c     based on computing the needed integrals via a reference quadrature.

         alf=alpha
         bet=beta
         r=1
         rwt=r+nq
         wtfn=rwt+nq
         scr=wtfn+nq
         ply=scr+nq
         arg=ply+nq*(nord+1)
         need=wpadti(arg+nq)
         call getmem(need,pxwt,ngot(1),'xwt',0)

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

         call vfill(xwt(ply),1.d0,nq)
         mu=scaprd(xwt(ply),xwt(ply),xwt(rwt),xwt(wtfn),nq)
         write(iout,7) mu
         if(parity.eq.'none') then 
            call copy(xwt(r),xwt(arg),nq)
         else
            if( dollar('$v0('//coord//')',card,title,inp) ) then
                angmom=intkey(card,'angular-momentum',0,' ')
            endif
            call vmul(xwt(arg),xwt(r),xwt(r),nq)
            ieo=angmom/2
            tst = angmom - 2*ieo
            if(tst.eq.1) then
               call copy(xwt(r),xwt(ply),nq)
            endif
            endpts(1)=endpts(1)*endpts(1)
            endpts(2)=endpts(2)*endpts(2)
         endif

c        Note that lancz will return the recursion coefficients on
c        the actual interval.  This is consistent with what is done for
c        the known cases.  The starting function is specified or taken
c        to be unity.

         call lancz(xwt(ply),xwt(arg),a,b,xwt(rwt),xwt(wtfn),
     1              xwt(scr),nq,nord)
  
c        Get rid of the unneeded memory

         call getmem(-ngot(1),pxwt,idum,'xwt',idum)

c        Modify them if using a Lobatto quadrature

         call modab(a,b,nfix,endpts,nord)
         if(prn(3)) then
            title='lanczos a coefficients'
            call prntfm(title,a,nord,1,nord,1,iout)
            title='lanczos b coefficients'
            call prntfm(title,b,n,1,n,1,iout)
         endif

c        Get the points and weights and then compute the coordinate
c        functions

         eigv=1
         scr=eigv+nord*nord
         need=wpadti(scr+npt)
         call getmem(need,py,ngot(3),'junk',0)
         call copy(a,pt,nord)
         call copy(b,y(scr),nord)       

c        Generate the non-classical points and weights and the
c        transformation matrix from the orthogonal polynomials
c        to the co-ordinate functions.

        call genq(pt,y(scr),wpt,y(eigv),fix,endpts,mu,nord)
         if(parity.ne.'none') then
            call xsq2x(pt,edge,nord)
         endif  
         title='eigenvectors'
         call prntrm(title,y(eigv),nord,nord,nord,nord,iout)
         if(prn(4)) then
            title='final nodes'
            call prntfm(title,pt,nord,1,nord,1,iout)
            title='final weights'
            call prntfm(title,wpt,nord,1,nord,1,iout)
         endif

c        Generate the needed functions at all required points.

         call cpoly(cp,dcp,ddcp,pt,y(scr),n,npt,parity,prn(5))
         call genwtf(pt,y(scr),wtt,dwtt,ddwtt,alpha,beta,
     1               edge,npt,typwt,prn(5))
         call getmem(-ngot(3),py,idum,'junk',idum)
      endif
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
 7    format(/,1x,'weight integral = ',e15.8)
      end
