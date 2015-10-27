*deck cfcoef.f
c***begin prologue     cfcoef
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           orthogonal polynomials
c***author             schneider, barry (nsf)
c***source             math
c***purpose            expand complex function in orthogonal polynomials.
c***                   
c***references         
c
c***routines called    
c***end prologue       cfcoef
      subroutine cfcoef(f,df,ddf,p,dp,ddp,x,wts,c,dum,n,npts)
      implicit integer (a-z)
      real*8 p, dp, ddp, x, wts
      complex*16 eye, f, df, ddf, c, dum
      character*80 title
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1) 
      dimension c(n), f(npts), df(npts), ddf(npts), x(npts), wts(npts)
      dimension dum(npts,3)
      common/io/inp, iout 
      data eye/(0.d0,1.d0)/
      do 10 i=1,npts
         f(i) = x(i) + eye*( exp(eye*x(i)) - 1.d0 )
         df(i) = 1.d0 - exp(eye*x(i))
         ddf(i)= -eye*exp(eye*x(i))
c          f(i)=x(i)*x(i)
c          df(i)=2.d0*x(i)
c          ddf(i)=2.d0 
 10   continue   
      do 20 i=1,npts
         f(i) = f(i)*wts(i)
 20   continue   
      call ebtcc(c,p,f,n,npts,1)
      title='polynomial coefficients'
      call prntcm(title,c,n,1,n,1,iout)
      do 30 i=1,npts
         f(i) = f(i)/wts(i)
 30   continue   
      call ebcc(dum(1,1),p,c,npts,n,1)      
      call ebcc(dum(1,2),dp,c,npts,n,1)      
      call ebcc(dum(1,3),ddp,c,npts,n,1)      
      do 40 i=1,npts
         write(iout,1) x(i)
         write(iout,2) f(i), dum(i,1)
         write(iout,3) df(i), dum(i,2)
         write(iout,4) ddf(i), dum(i,3)
 40   continue   
      return
 1    format(/,'          x = ',e15.8)
 2    format('fexact    = ','(',e15.8,',',e15.8,')'/,
     1       'fapprox   = ','(',e15.8,',',e15.8,')')
 3    format('dfexact   = ','(',e15.8,',',e15.8,')'/,
     1       'dfapprox  = ','(',e15.8,',',e15.8,')')
 4    format('ddfexact  = ','(',e15.8,',',e15.8,')'/,
     1       'ddfapprox = ','(',e15.8,',',e15.8,')')
      end       
