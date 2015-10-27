*deck fcoef.f
c***begin prologue     fcoef
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           orthogonal polynomials
c***author             schneider, barry (nsf)
c***source             math
c***purpose            expand function in orthogonal polynomials.
c***                   
c***references         
c
c***routines called    
c***end prologue       fcoef
      subroutine fcoef(f,df,ddf,p,dp,ddp,x,wts,c,dum,n,npts,nocoef)
      implicit integer (a-z)
      real*8 p, dp, ddp, c, f, df, ddf, x, wts, dum
      character*80 title
      logical nocoef
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1) 
      dimension c(n), f(npts), df(npts), ddf(npts), x(npts), wts(npts)
      dimension dum(npts,3)
      common/io/inp, iout 
      do 10 i=1,npts
c         f(i)=x(i)*exp(-5.d0*x(i))
c         df(i)=-5.d0*x(i)*exp(-5.d0*x(i))+exp(-5.d0*x(i))
c         ddf(i)=25.d0*x(i)*exp(-5.d0*x(i))-10.d0*exp(-5.d0*x(i))
c         f(i)=sin(5.d0*x(i))
c         df(i)=5.d0*cos(5.d0*x(i))
c         ddf(i)=-25.d0*sin(5.d0*x(i))
c          f(i)=sqrt(x(i))
c          df(i)=.5d0/f(i)
c          ddf(i)=-.25d0/(x(i)*f(i))
          f(i)=x(i)*x(i)*x(i)
          df(i)=3.d0*x(i)*x(i)
          ddf(i)=6.d0*x(i)
 10   continue   
      if(.not.nocoef) then
         call vmul(f,f,wts,npts)
         call ebtc(c,p,f,n,npts,1)
         title='polynomial coefficients'
         call prntrm(title,c,n,1,n,1,iout)
         call vdiv(f,f,wts,npts)
      endif
      call ebc(dum(1,1),p,c,npts,n,1)      
      call ebc(dum(1,2),dp,c,npts,n,1)      
      call ebc(dum(1,3),ddp,c,npts,n,1)      
      do 20 i=1,npts
         write(iout,1) x(i)
         write(iout,2) f(i), dum(i,1)
         write(iout,3) df(i), dum(i,2)
         write(iout,4) ddf(i), dum(i,3)
 20   continue   
      return
 1    format(/,'          x = ',e15.8)
 2    format('fexact   = ',e15.8,1x,'fapprox   = ',e15.8)
 3    format('dfexact  = ',e15.8,1x,'dfapprox  = ',e15.8)
 4    format('ddfexact = ',e15.8,1x,'ddfapprox = ',e15.8)
      end       
