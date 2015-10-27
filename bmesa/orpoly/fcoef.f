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
      subroutine fcoef(f,df,ddf,pn,dpn,ddpn,x,wts,c,dum,n,npts,nocoef)
      implicit integer (a-z)
      real*8 pn, dpn, ddpn, c, f, df, ddf, x, wts, dum
      character*80 title
      logical nocoef
      dimension pn(npts,0:n-1), dpn(npts,0:n-1), ddpn(npts,0:n-1) 
      dimension c(n), f(npts), df(npts), ddf(npts), x(npts), wts(n)
      dimension dum(npts,3)
      common/io/inp, iout 
      do 10 i=1,npts
c         f(i)=sin(x(i))
c         df(i)=cos(x(i))
c         ddf(i)=-f(i)
          f(i)=x(i)
          df(i)=1.d0
          ddf(i)=0.d0
 10   continue   
      if(.not.nocoef) then
         call vmul(f,f,wts,n)
         call ebtc(c,pn,f,n,n,1)
         title='polynomial coefficients'
         call prntrm(title,c,n,1,n,1,iout)
         call vdiv(f,f,wts,n)
      endif
      call ebc(dum(1,1),pn,c,npts,n,1)      
      call ebc(dum(1,2),dpn,c,npts,n,1)      
      call ebc(dum(1,3),ddpn,c,npts,n,1)      
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
