c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Generalized Orthogonal Polynomials and Derivatives}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle    

*deck poly.f
c***begin prologue     poly
c***date written       022202   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           coordinate eigenfunctions
c***author             schneider, barry (nsf)
c***source
c***purpose            generate orthogonal polynomials using recursion.
c***                   compute coordinate functions by transformation
c***                   to representation diagonalizing coordinate.
c***           
c***description
c***references         
c
c***routines called    
c***end prologue       poly

c\begin{eqnarray}
c  \beta_{j} P_{j}(x) &=& \big ( x - \alpha_{j} \big ) P_{j-1}(x) - \beta_{j-1} P_{j-2}(x) 
c                             \\ \nonumber
c  \alpha_{j} &=& \langle P_{j-1} \mid x \mid P_{j-1} \rangle \\ \nonumber
c  \beta_{j} &=& \langle P_{j-1} \mid x \mid P_{j} \rangle \\ \nonumber
c  \beta_{j} P_{j}^{\prime}(x) &=& \big ( x - \alpha_{j} \big ) P_{j-1}^{\prime}(x) 
c                             - \beta_{j-1} P_{j-2}^{\prime}(x) + P_{j-1}(x)
c                                 \\ \nonumber
c  \beta_{j} P_{j}^{\prime \prime}(x) &=& \big ( x - \alpha_{j} \big ) 
c                                      P_{j-1}^{\prime \prime}(x) 
c                            - \beta_{j-1} P_{j-2}^{\prime \prime}(x) 
c                            + 2 P_{j-1}^{\prime}(x)
c\end{eqnarray} 

      subroutine poly(p,dp,ddp,cp,dcp,ddcp,pt,arg,a,b,eigv,
     1                mu,n,npt,typ,prn)
      implicit integer (a-z)
      real*8 p, dp, ddp, cp, dcp, ddcp, pt, a, b, eigv, mu
      real*8 arg, fac, fac2, tmp
      character*(*) typ
      character*80 title
      logical prn
      dimension p(npt,0:n), dp(npt,0:n), ddp(npt,0:n)
      dimension cp(npt,0:n), dcp(npt,0:n), ddcp(npt,0:n)
      dimension pt(npt), arg(npt), a(1:n), b(1:n), eigv(0:n,0:n)
      dimension prn(2)
      common/io/inp, iout 

c     Using the analytic or computed recursion coefficients, $\alpha$ and
c     $\beta$, we compute the orthogonal polynomials and the coordinate
c     functions for arbitrary points.

      if(typ.eq.'quadratic') then
        call vmul(arg,pt,pt,npt)
      else
        call copy(pt,arg,npt)
      endif

c     Do the first and second functions as special cases.

      fac=1.d0/sqrt(mu)
      call vfill(p(1,0),fac,npt)
      do 10 i=1,npt
         p(i,1) = (arg(i) - a(1))*p(i,0)
	 p(i,1) = p(i,1)/b(1)
 10   continue
      
c     Recur for the remaining ones

      do 20 j=2,n
         do 30 i=1,npt
            p(i,j) = (arg(i) - a(j) )*p(i,j-1) - b(j-1)*p(i,j-2)
            p(i,j) = p(i,j)/b(j)
 30      continue
 20   continue   

c     Lets do the derivatives

      call rzero(dp(1,0),npt)
      do 40 i=1,npt
         dp(i,1) = p(i,0)/b(1)
 40   continue
      do 50 j=2,n
         do 60 i=1,npt
            dp(i,j) = ( arg(i) - a(j) )*dp(i,j-1) - b(j-1)*dp(i,j-2)
     1                     + p(i,j-1)
            dp (i,j) = dp(i,j)/b(j)
 60      continue
 50   continue
      call rzero(ddp(1,0),npt)   
      call rzero(ddp(1,1),npt)   
      do 70 j=2,n
         do 80 i=1,npt
            ddp(i,j) = ( arg(i) - a(j) )*ddp(i,j-1) - b(j-1)*ddp(i,j-2)
     1                    + 2.d0*dp(i,j-1)
            ddp(i,j) = ddp(i,j)/b(j)
 80      continue
 70   continue
c     These are derivatives wrt arg.  Convert to derivatives wrt pt
c      if(typ.eq.'quadratic') then
c         do 100 j=1,npt
c            fac=2.d0*pt(j)
c            fac2=fac*fac
c            do 200 i=0,n
c               tmp=dp(j,i)
c               dp(j,i) = fac*tmp
c               ddp(j,i) = 2.d0*tmp + fac2*ddp(j,i)
c 200        continue
c 100     continue   
c      endif
      if(prn(1)) then
         title='orthogonal polynomials'
         call prntfm(title,p(1,0),npt,n+1,npt,n+1,iout)   
         title='first derivative of orthogonal polynomials'
         call prntfm(title,dp(1,0),npt,n+1,npt,n+1,iout)   
         title='second derivative of orthogonal polynomials'
         call prntfm(title,ddp(1,0),npt,n+1,npt,n+1,iout)   
      endif

c     Transform to the co-ordinate functions
     
      call ebc(cp,p,eigv,npt,n+1,n+1)
      call ebc(dcp,dp,eigv,npt,n+1,n+1)
      call ebc(ddcp,ddp,eigv,npt,n+1,n+1)
      if(prn(2)) then
         title='coordinate function'
         call prntfm(title,cp(1,0),npt,n+1,npt,n+1,iout)   
         title='first derivative of coordinate function'
         call prntfm(title,dcp(1,0),npt,n+1,npt,n+1,iout)   
         title='second derivative of coordinate function'
         call prntfm(title,ddcp(1,0),npt,n+1,npt,n+1,iout)   
      endif 
      return
      end       
