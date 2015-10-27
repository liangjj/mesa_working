*deck @(#)qpow.f	2.1  10/10/91
      function qpow (n,l,alpha,rk,t,expab)
      implicit real*8(a-h,o-z)
c     power series for q(n,l).......
c
c     ----- functions returned -----
      real*8 qpow
c     ----- arguments unchanged -----
      integer n,l
      real*8 alpha,rk,t,expab
c
c     ----- local variables -----
      real*8 zero,one,two,tenm14,round
      real*8 xkp,prefac,term,sum,xj,xnum,xden
c     ----- common -----
      real*8 dfac,pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/dfac/dfac(23)
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00)
      parameter (tenm14=1.0d-14,round=0.4d-14)
c
c
      if (l.eq.0) then
         xkp=one
      else
         xkp=rk**l
      endif
      prefac=expab*xkp/sqrt((two*alpha)**(n+l+1))
      xnum=l+n-1
      xden=l+l+1
      term=dfac(l+n+1)/dfac(l+l+3)
      sum=term
      xj=zero
   10 xnum=xnum+two
         xden=xden+two
         xj=xj+one
         term=term*t*xnum/(xj*xden)
         sum=sum+term
      if ((term/sum).gt.tenm14) go to 10
c
      qpow=prefac*sum
      if (mod((l+n),2).eq.0) qpow=qpow*sqpi2
c     following stmt corrects for truncation error.
      qpow=qpow*(one+xj*round)
c
c
      return
      end
